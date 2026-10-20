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

file_name = "snail_mucin"

mindt = 0.3 # Time interval between consecutive frames in seconds
pxsz = 0.1465 # Pixel size in micrometers per pixel

# Load the microscopy video
intensity=suppressWarnings(readTIF(paste0("../data/", file_name, ".tif")))
# Estimate the MSD, its 95% confidence interval, and the storage and loss moduli
M = 200 # Number of tracer particles
sam_mf = SAM(intensity = intensity, intensity_str = "SST_array", M = M, mindt = mindt, pxsz = pxsz,
             model_free_estimation = TRUE, uncertainty = TRUE, estimate_moduli = TRUE)

# MF-AIUQ results
d_input = sam_mf@d_input
MSD_MFAIUQ = sam_mf@msd_est
MSD_lower = sam_mf@msd_lower
MSD_upper = sam_mf@msd_upper
# Storage and Loss moduli
GSER_MF_AIUQ = data.frame(omega = sam_mf@omega, Gp = sam_mf@Gp, Gpp = sam_mf@Gpp)


# Plot
par(mar = c(3.9, 5.6, 1.5, 1.5), mgp = c(3, 1, 0), xpd = NA)
xlim_l = min(d_input[-1])
xlim_u = max(d_input[-1]) + 2
ylim_l = min(0.001, min(MSD_lower[-1])) * 0.4
ylim_u = max(range((c(MSD_MFAIUQ[-1])))[2], max(MSD_upper[-1])) * 1.2
plot(d_input[-1], MSD_MFAIUQ[-1], log = "xy", type = "n",
     xlab = expression(Delta*t ~ "[s]"), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">" ~ "[" * mu*m^2 * "]"), ,
     xlim = c(xlim_l, xlim_u),
     ylim = c(ylim_l, ylim_u),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", cex.lab = 2.2)
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
points(d_input[-1], MSD_MFAIUQ[-1], pch = 1, col = "cyan", cex = 1.8, lwd = 1.8)
polygon(c(d_input[-1], rev(d_input[-1])),
        c(MSD_lower[-1], rev(MSD_upper[-1])),
        col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)


# G' and G'' MF-AIUQ
par(mar = c(3.9, 5.6, 1.5, 1.5), mgp = c(3, 1, 0))
xlim_l = 0.001
xlim_u = max(GSER_MF_AIUQ$omega) * 2
ylim_l = 0.001
ylim_u = max(GSER_MF_AIUQ$Gp, GSER_MF_AIUQ$Gpp) * 2
plot(GSER_MF_AIUQ$omega, GSER_MF_AIUQ$Gp, log = "xy", type = "l", col = "cyan", lwd = 5, 
     xlab = expression(omega ~ "[" * rad/s * "]"), ylab = expression(G*minute ~ "," ~ G*minute*minute ~ "[Pa]"),
     xlim = c(xlim_l, xlim_u), ylim = c(ylim_l, ylim_u),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", cex.lab = 2)
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
points(GSER_MF_AIUQ$omega, GSER_MF_AIUQ$Gpp,
       pch = 21, bg = "white", col = "cyan", cex = 1.3)
legend("bottomright",
       legend = c(expression(MF-AIUQ~G*minute),
                  expression(MF-AIUQ~G*minute*minute)),
       col = c("cyan", "cyan"),
       lty = c(1, NA),
       pch = c(NA, 21),
       lwd = c(2, 1.5),
       pt.bg = c(NA, "white"),
       pt.cex = 1.3,
       bty = "n",
       cex = 1.2)

