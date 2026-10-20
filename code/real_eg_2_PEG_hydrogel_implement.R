# Load required packages
library(AIUQ)
library(bioimagetools)

# The working directory should be the code/ folder in the downloaded MF-AIUQ repository.
# Check the current working directory with:
# getwd()

# If necessary, replace the example path below with the location of the
# code/ folder on your computer.

# macOS or Linux:
# setwd("/path/to/MF-AIUQ/code")

# Windows:
# setwd("C:/path/to/MF-AIUQ/code")

# Expected directory structure:
# MF-AIUQ/
# ├── code/    <- working directory
# ├── data/
# └── results/

# Download the microscopy videos from Zenodo (link) and place the TIFF files
# in the data/ folder.

# Select one scenario:
# "13min", "22min", "31min", "40min", "49min", or "58min"
scenario = "13min"

# Load the microscopy video
intensity = suppressWarnings(readTIF(paste0("../data/20mgmL_25C_DMSO_", scenario, ".tif")))

mindt = 0.0309 # Time interval between consecutive frames in seconds
pxsz = 0.2930 # Pixel size in micrometers
M = 200 # Number of tracer particles

# Estimate the MSD and its 95% confidence interval using MF-AIUQ
sam_mf = SAM(intensity = intensity, intensity_str = "SST_array", M = M, mindt = mindt, pxsz = pxsz,
             model_free_estimation = TRUE, uncertainty = TRUE)

# MF-AIUQ results
d_input = sam_mf@d_input
MSD_MFAIUQ = sam_mf@msd_est
MSD_lower = sam_mf@msd_lower
MSD_upper = sam_mf@msd_upper


# Plot
if (scenario == "13min") {
  legend_text = "13 min"
} else if (scenario == "22min") {
  legend_text = "22 min"
} else if (scenario == "31min") {
  legend_text = "31 min"
} else if (scenario == "40min") {
  legend_text = "40 min"
} else if (scenario == "49min") {
  legend_text = "49 min"
} else if (scenario == "58min") {
  legend_text = "58 min"
}

par(mar = c(4.1, 5.6, 1.5, 1.5), mgp = c(3.2, 1, 0), xpd = NA)
xlim_l = 0.0309
xlim_u = max(d_input[-1]) + 2
ylim_l = min(0.001, min(MSD_lower[-1])) * 0.9
ylim_u = max(range((c(MSD_MFAIUQ[-1])))[2], max(MSD_upper[-1])) * 1.2
plot(d_input[-1], MSD_MFAIUQ[-1], log = "xy", pch = 1, col = "cyan",
     xlab = expression(Delta*t ~ "[s]"), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">" ~ "[" * mu*m^2 * "]"),
     xlim = c(xlim_l, xlim_u),
     ylim = c(ylim_l, ylim_u),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", cex.lab = 1.8)
kx = floor(log10(c(xlim_l, xlim_u)))
major_x = 10^(kx[1]:kx[2])
major_x = major_x[major_x >= xlim_l & major_x <= xlim_u]
minor_x = as.vector(outer(10^(kx[1]:kx[2]), 2:9))
minor_x = minor_x[minor_x >= xlim_l & minor_x <= xlim_u]
axis(1, at = major_x, labels = parse(text = paste0("10^", log10(major_x))), tcl = 0.4, lwd = 1.4, cex.axis = 1.8)
axis(1, at = minor_x, labels = FALSE, tcl = 0.3, lwd = 1)
ky = floor(log10(c(ylim_l, ylim_u)))
major_y = 10^(ky[1]:ky[2])
major_y = major_y[major_y >= ylim_l & major_y <= ylim_u]
minor_y = as.vector(outer(10^(ky[1]:ky[2]), 2:9))
minor_y = minor_y[minor_y >= ylim_l & minor_y <= ylim_u]
axis(2, at = major_y, labels = parse(text = paste0("10^", log10(major_y))), tcl = 0.4, lwd = 1.4, cex.axis = 1.8)
axis(2, at = minor_y, labels = FALSE, tcl = 0.3, lwd = 1)
polygon(c(d_input[-1], rev(d_input[-1])),
        c(MSD_lower[-1], rev(MSD_upper[-1])),
        col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)
legend("bottomright", legend = legend_text,
       bty = "n", text.font = 1.8, cex = 2.2, inset = c(-0.03, -0.04))
