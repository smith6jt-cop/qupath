/**
 * Radial Line CD20 Intensity Profiler for Follicle Identification
 *
 * Draws lines through a region at 5-degree increments (0° to 175°),
 * measures DAB (CD20) intensity along each line, and identifies
 * regions of high CD20 expression as potential follicles.
 *
 * Usage:
 *   1. Open an H-DAB stained image in QuPath with CD20 marker
 *   2. Optionally select an annotation to limit the analysis area
 *   3. Run this script in the QuPath Script Editor (Automate > Script Editor)
 *   4. Results: Follicle annotations are created in the object hierarchy
 *
 * How it works:
 *   - Lines radiate from the center of the image (or selected annotation)
 *   - Starting vertical (0°), rotating 5° clockwise each step through 175°
 *   - DAB optical density is sampled along each line using color deconvolution
 *   - Contiguous segments above the DAB threshold are identified as CD20+ regions
 *   - CD20+ points from all lines are spatially clustered into follicle candidates
 *   - An elliptical annotation is created around each cluster
 *
 * Requirements:
 *   - Image must have H-DAB color deconvolution stains set
 *     (Image > Set image type > Brightfield (H-DAB))
 */

import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.RegionRequest
import qupath.lib.color.ColorTransformer
import qupath.lib.color.ColorDeconvolutionHelper
import java.awt.image.BufferedImage

// ====================== CONFIGURATION ======================
// Rotation parameters
double ANGLE_STEP_DEG   = 5       // Degrees between each radial line
double START_ANGLE_DEG  = 0       // Starting angle (0 = vertical)
double END_ANGLE_DEG    = 175     // End angle (exclusive of 180, which duplicates 0)

// Intensity sampling
double SAMPLE_SPACING   = 2.0     // Pixels between intensity samples along each line
double DOWNSAMPLE       = 1.0     // Downsample factor for reading pixels (increase for speed on large images)

// CD20+ detection thresholds
double DAB_THRESHOLD    = 0.25    // Minimum DAB optical density to consider CD20+
int    MIN_SEGMENT_PX   = 50      // Minimum length (in pixels) of continuous CD20+ signal
int    SMOOTHING_WINDOW = 11      // Sliding window size for smoothing intensity profiles (odd number)

// Follicle clustering
double CLUSTER_RADIUS   = 150.0   // Max distance (px) between points in the same follicle cluster
int    MIN_CLUSTER_PTS  = 3       // Minimum number of segment midpoints to form a follicle
double FOLLICLE_PADDING = 1.3     // Multiplier to pad the ellipse around each cluster

// Output options
boolean ADD_SCAN_LINES       = false  // Add the radial scan lines as annotations (for visualization)
boolean ADD_FOLLICLE_REGIONS = true   // Add detected follicle ellipses as annotations
String  FOLLICLE_CLASS_NAME  = "Follicle"  // PathClass name for follicle annotations
// ============================================================


// ---------- Helper: Smooth an array with a sliding mean ----------
double[] smooth(double[] values, int window) {
    if (window < 2) return values
    int halfW = (int)(window / 2)
    double[] result = new double[values.length]
    for (int i = 0; i < values.length; i++) {
        int lo = Math.max(0, i - halfW)
        int hi = Math.min(values.length - 1, i + halfW)
        double sum = 0
        for (int j = lo; j <= hi; j++) sum += values[j]
        result[i] = sum / (hi - lo + 1)
    }
    return result
}

// ---------- Helper: Simple DBSCAN-style spatial clustering ----------
List<List<double[]>> clusterPoints(List<double[]> points, double eps, int minPts) {
    int n = points.size()
    int[] labels = new int[n]
    Arrays.fill(labels, -1)  // -1 = unvisited
    int clusterId = 0

    for (int i = 0; i < n; i++) {
        if (labels[i] != -1) continue
        // Find neighbors
        List<Integer> neighbors = []
        for (int j = 0; j < n; j++) {
            double dx = points[i][0] - points[j][0]
            double dy = points[i][1] - points[j][1]
            if (Math.sqrt(dx * dx + dy * dy) <= eps)
                neighbors.add(j)
        }
        if (neighbors.size() < minPts) {
            labels[i] = -2  // noise
            continue
        }
        labels[i] = clusterId
        List<Integer> seedSet = new ArrayList<>(neighbors)
        seedSet.remove(Integer.valueOf(i))
        int idx = 0
        while (idx < seedSet.size()) {
            int q = seedSet[idx]
            if (labels[q] == -2) labels[q] = clusterId  // noise becomes border
            if (labels[q] != -1) { idx++; continue }
            labels[q] = clusterId
            List<Integer> qNeighbors = []
            for (int j = 0; j < n; j++) {
                double dx = points[q][0] - points[j][0]
                double dy = points[q][1] - points[j][1]
                if (Math.sqrt(dx * dx + dy * dy) <= eps)
                    qNeighbors.add(j)
            }
            if (qNeighbors.size() >= minPts) {
                for (int nb : qNeighbors) {
                    if (!seedSet.contains(nb)) seedSet.add(nb)
                }
            }
            idx++
        }
        clusterId++
    }

    // Group points by cluster
    Map<Integer, List<double[]>> clusters = [:]
    for (int i = 0; i < n; i++) {
        if (labels[i] >= 0) {
            clusters.computeIfAbsent(labels[i], { [] }).add(points[i])
        }
    }
    return clusters.values().toList()
}


// ==================== MAIN SCRIPT ====================
def imageData = getCurrentImageData()
if (imageData == null) {
    println "ERROR: No image is open. Please open an H-DAB image first."
    return
}

def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()
def stains = imageData.getColorDeconvolutionStains()

if (stains == null) {
    println "ERROR: No color deconvolution stains are set."
    println "Please set the image type to 'Brightfield (H-DAB)' via Image > Set image type."
    return
}

println "Stain 1: ${stains.getStain(1).getName()}"
println "Stain 2: ${stains.getStain(2).getName()}"

// Determine analysis region (selected annotation or full image)
def selectedObject = getSelectedObject()
double regionX, regionY, regionW, regionH
if (selectedObject != null && selectedObject.getROI() != null) {
    def roi = selectedObject.getROI()
    regionX = roi.getBoundsX()
    regionY = roi.getBoundsY()
    regionW = roi.getBoundsWidth()
    regionH = roi.getBoundsHeight()
    println "Using selected annotation as region of interest"
} else {
    regionX = 0
    regionY = 0
    regionW = server.getWidth()
    regionH = server.getHeight()
    println "No selection found - analyzing entire image"
}

double centerX = regionX + regionW / 2.0
double centerY = regionY + regionH / 2.0
double maxRadius = Math.sqrt(regionW * regionW + regionH * regionH) / 2.0

println "Region: (${regionX}, ${regionY}) size ${regionW} x ${regionH}"
println "Center: (${centerX}, ${centerY}), Max radius: ${maxRadius}"

// Precompute OD lookup tables for faster pixel processing
double[] odLutR = ColorDeconvolutionHelper.makeODLUT(stains.getMaxRed())
double[] odLutG = ColorDeconvolutionHelper.makeODLUT(stains.getMaxGreen())
double[] odLutB = ColorDeconvolutionHelper.makeODLUT(stains.getMaxBlue())
double[][] invMat = stains.getMatrixInverse()

// DAB is stain index 1 (column 1 in the inverse matrix)
int dabCol = 1

// Collect all CD20+ segment midpoints across all lines
List<double[]> allSegmentMidpoints = []
List<Map> allSegments = []  // store full segment info for reporting

int nLines = (int)((END_ANGLE_DEG - START_ANGLE_DEG) / ANGLE_STEP_DEG) + 1
println "\nScanning ${nLines} radial lines (${START_ANGLE_DEG}° to ${END_ANGLE_DEG}° in ${ANGLE_STEP_DEG}° steps)..."

def newAnnotations = []

for (int li = 0; li < nLines; li++) {
    double angleDeg = START_ANGLE_DEG + li * ANGLE_STEP_DEG
    double angleRad = Math.toRadians(angleDeg)

    // Line direction: 0° = vertical (pointing up), rotating clockwise
    // So dx = sin(angle), dy = -cos(angle)
    double dx = Math.sin(angleRad)
    double dy = -Math.cos(angleRad)

    // Line endpoints extending from center to region edges
    double x1 = centerX - dx * maxRadius
    double y1 = centerY - dy * maxRadius
    double x2 = centerX + dx * maxRadius
    double y2 = centerY + dy * maxRadius

    // Clip endpoints to image bounds
    x1 = Math.max(0, Math.min(server.getWidth() - 1, x1))
    y1 = Math.max(0, Math.min(server.getHeight() - 1, y1))
    x2 = Math.max(0, Math.min(server.getWidth() - 1, x2))
    y2 = Math.max(0, Math.min(server.getHeight() - 1, y2))

    // Optionally add scan line annotation
    if (ADD_SCAN_LINES) {
        def lineROI = ROIs.createLineROI(x1, y1, x2, y2, server.getMetadata().getDefaultPlane())
        def lineAnnotation = PathObjects.createAnnotationObject(lineROI, PathClass.fromString("Scan Line"))
        newAnnotations.add(lineAnnotation)
    }

    // Sample DAB intensity along the line
    double lineLength = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
    int nSamples = (int)(lineLength / SAMPLE_SPACING)
    if (nSamples < 2) continue

    double stepX = (x2 - x1) / nSamples
    double stepY = (y2 - y1) / nSamples

    double[] dabProfile = new double[nSamples]
    double[] sampleXcoords = new double[nSamples]
    double[] sampleYcoords = new double[nSamples]

    // Read pixels in small tiles along the line
    for (int s = 0; s < nSamples; s++) {
        double sx = x1 + s * stepX
        double sy = y1 + s * stepY
        sampleXcoords[s] = sx
        sampleYcoords[s] = sy

        int px = (int) Math.round(sx)
        int py = (int) Math.round(sy)

        // Bounds check
        if (px < 0 || px >= server.getWidth() || py < 0 || py >= server.getHeight()) {
            dabProfile[s] = 0
            continue
        }

        // Read a 1x1 pixel region
        try {
            def request = RegionRequest.createInstance(server.getPath(), DOWNSAMPLE, px, py, 1, 1)
            BufferedImage img = server.readRegion(request)
            int rgb = img.getRGB(0, 0)

            // Color deconvolution to get DAB optical density
            double rOD = odLutR[(rgb >> 16) & 0xFF]
            double gOD = odLutG[(rgb >> 8) & 0xFF]
            double bOD = odLutB[rgb & 0xFF]
            dabProfile[s] = rOD * invMat[0][dabCol] + gOD * invMat[1][dabCol] + bOD * invMat[2][dabCol]
        } catch (Exception e) {
            dabProfile[s] = 0
        }
    }

    // Smooth the profile to reduce noise
    double[] smoothed = smooth(dabProfile, SMOOTHING_WINDOW)

    // Find contiguous segments above threshold
    boolean inSegment = false
    int segStart = 0
    for (int s = 0; s <= nSamples; s++) {
        double val = (s < nSamples) ? smoothed[s] : 0
        if (!inSegment && val >= DAB_THRESHOLD) {
            inSegment = true
            segStart = s
        } else if (inSegment && val < DAB_THRESHOLD) {
            inSegment = false
            int segEnd = s - 1
            int segLengthPx = (int)((segEnd - segStart) * SAMPLE_SPACING)

            if (segLengthPx >= MIN_SEGMENT_PX) {
                int mid = (segStart + segEnd) / 2
                double midX = sampleXcoords[mid]
                double midY = sampleYcoords[mid]

                // Compute mean DAB OD for this segment
                double segSum = 0
                for (int k = segStart; k <= segEnd; k++) segSum += smoothed[k]
                double segMean = segSum / (segEnd - segStart + 1)

                allSegmentMidpoints.add([midX, midY] as double[])
                allSegments.add([
                    angle: angleDeg,
                    startX: sampleXcoords[segStart], startY: sampleYcoords[segStart],
                    endX: sampleXcoords[segEnd], endY: sampleYcoords[segEnd],
                    midX: midX, midY: midY,
                    lengthPx: segLengthPx,
                    meanDAB: segMean
                ])
            }
        }
    }

    if ((li + 1) % 10 == 0 || li == nLines - 1) {
        println "  Processed line ${li + 1}/${nLines} (${angleDeg}°)"
    }
}

println "\nFound ${allSegments.size()} CD20+ segments across all lines"

// Print segment details
if (allSegments.size() > 0) {
    println "\n--- CD20+ Segment Summary ---"
    println String.format("%-8s %-10s %-12s %-12s", "Angle", "Length(px)", "Mean DAB OD", "Midpoint")
    for (seg in allSegments) {
        println String.format("%-8.1f %-10d %-12.3f (%.0f, %.0f)",
            seg.angle, seg.lengthPx, seg.meanDAB, seg.midX, seg.midY)
    }
}

// Cluster segment midpoints to identify follicles
if (allSegmentMidpoints.size() >= MIN_CLUSTER_PTS && ADD_FOLLICLE_REGIONS) {
    println "\nClustering ${allSegmentMidpoints.size()} segment midpoints (radius=${CLUSTER_RADIUS}, minPts=${MIN_CLUSTER_PTS})..."

    def clusters = clusterPoints(allSegmentMidpoints, CLUSTER_RADIUS, MIN_CLUSTER_PTS)
    println "Identified ${clusters.size()} follicle candidate(s)"

    def follicleClass = PathClass.fromString(FOLLICLE_CLASS_NAME)
    int follicleNum = 0

    for (cluster in clusters) {
        follicleNum++

        // Compute bounding ellipse for the cluster
        double sumX = 0, sumY = 0
        double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE
        double minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE
        for (pt in cluster) {
            sumX += pt[0]; sumY += pt[1]
            minX = Math.min(minX, pt[0]); maxX = Math.max(maxX, pt[0])
            minY = Math.min(minY, pt[1]); maxY = Math.max(maxY, pt[1])
        }
        double cx = sumX / cluster.size()
        double cy = sumY / cluster.size()
        double rx = (maxX - minX) / 2.0 * FOLLICLE_PADDING
        double ry = (maxY - minY) / 2.0 * FOLLICLE_PADDING

        // Ensure a minimum size
        rx = Math.max(rx, MIN_SEGMENT_PX)
        ry = Math.max(ry, MIN_SEGMENT_PX)

        def ellipseROI = ROIs.createEllipseROI(
            cx - rx, cy - ry, rx * 2, ry * 2,
            server.getMetadata().getDefaultPlane()
        )
        def follicleAnnotation = PathObjects.createAnnotationObject(ellipseROI, follicleClass)
        follicleAnnotation.setName("Follicle ${follicleNum}")

        // Add measurements to the annotation
        follicleAnnotation.measurements.put("Cluster size", (double) cluster.size())
        follicleAnnotation.measurements.put("Center X", cx)
        follicleAnnotation.measurements.put("Center Y", cy)

        // Compute mean DAB across all segments in this cluster
        double clusterDABsum = 0
        int clusterSegCount = 0
        for (pt in cluster) {
            for (seg in allSegments) {
                if (Math.abs(seg.midX - pt[0]) < 1 && Math.abs(seg.midY - pt[1]) < 1) {
                    clusterDABsum += seg.meanDAB
                    clusterSegCount++
                    break
                }
            }
        }
        if (clusterSegCount > 0) {
            follicleAnnotation.measurements.put("Mean DAB OD", clusterDABsum / clusterSegCount)
        }

        newAnnotations.add(follicleAnnotation)
        println "  Follicle ${follicleNum}: center=(${String.format('%.0f', cx)}, ${String.format('%.0f', cy)}), " +
                "radius=(${String.format('%.0f', rx)} x ${String.format('%.0f', ry)}), " +
                "${cluster.size()} line crossings"
    }
} else if (allSegmentMidpoints.size() < MIN_CLUSTER_PTS) {
    println "\nToo few CD20+ segments found (${allSegmentMidpoints.size()}) to identify follicles."
    println "Try adjusting DAB_THRESHOLD (currently ${DAB_THRESHOLD}) or MIN_SEGMENT_PX (currently ${MIN_SEGMENT_PX})."
}

// Add all new annotations to the hierarchy
if (!newAnnotations.isEmpty()) {
    hierarchy.addObjects(newAnnotations)
    hierarchy.getSelectionModel().clearSelection()
    println "\nAdded ${newAnnotations.size()} annotation(s) to the image"
}

println "\n=== Analysis Complete ==="
println "Tip: If no follicles were detected, try lowering DAB_THRESHOLD or MIN_SEGMENT_PX."
println "Tip: Set ADD_SCAN_LINES = true to visualize the radial scan pattern."
println "Tip: Adjust CLUSTER_RADIUS if follicles are being merged or split incorrectly."
