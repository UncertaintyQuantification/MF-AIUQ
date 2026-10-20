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

# Select one experiment:
# "1 um", "200 nm", "100 nm", or "100 nm BF"
scenario = "1 um"

if (scenario == "1 um") {
  file_name = "4pct_PVA_1um_25C"
} else if (scenario == "200 nm") {
  file_name = "4pct_PVA_200nm_25C"
} else if (scenario == "100 nm") {
  file_name = "4pct_PVA_100nm_25C"
} else if (scenario == "100 nm BF") {
  file_name = "4pct_PVA_100nm_25C_BF"
}
# Load the microscopy video from the data/ folder
intensity=suppressWarnings(readTIF(paste0("../data/", file_name, ".tif")))

mindt = 0.0309 # Time interval between consecutive frames in seconds
pxsz =0.2930 # Pixel size in micrometers

# Estimate the MSD and its 95% confidence interval using MF-AIUQ
M = 100 # Number of tracer particles
sam_mf = SAM(intensity = intensity, intensity_str = "SST_array", M = M, mindt = mindt, pxsz = pxsz,
             model_free_estimation = TRUE, uncertainty = TRUE)

# MF-AIUQ results
d_input = sam_mf@d_input
MSD_MFAIUQ = sam_mf@msd_est
MSD_lower = sam_mf@msd_lower
MSD_upper = sam_mf@msd_upper

if (scenario == "1 um") {
  file_name = "4pct_PVA_1um_25C"
  b_truth = 0.0759
  plot_name = "1um"
  ylim_l = min(MSD_lower[-1]) * 0.7
  legend_text = latex2exp::TeX("1 $\\mu$m")
} else if (scenario == "200 nm") {
  file_name = "4pct_PVA_200nm_25C"
  b_truth = 0.3795
  plot_name = "200nm"
  ylim_l = min(MSD_lower[-1]) * 0.1
  legend_text = "200 nm"
} else if (scenario == "100 nm") {
  file_name = "4pct_PVA_100nm_25C"
  b_truth = 0.7589
  plot_name = "100nm"
  ylim_l = min(MSD_lower[-1]) * 0.05
  legend_text = "100 nm"
} else if (scenario == "100 nm BF") {
  file_name = "4pct_PVA_100nm_25C_BF"
  b_truth = 0.7589
  plot_name = "100nm_BF"
  ylim_l = min(MSD_lower[-1]) * 0.05
  legend_text = "100 nm BF"
}

# Reference MSD
truth = d_input*b_truth


# Plot
img = intensity[,,250]
img = (img-min(img))/diff(range(img))
par(mar = c(4.1, 5.6, 1.5, 1.5), mgp = c(3.2, 1, 0), xpd = NA)
xlim_l = 0.0309
xlim_u = max(d_input[-1]) + 2
ylim_u = max(range((c(truth[-1], MSD_MFAIUQ[-1])))[2], max(MSD_upper[-1])) * 1.1
plot(d_input[-1], truth[-1], log = "xy", type = "n",
     xlab = expression(Delta*t ~ "[s]"), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">" ~ "[" * mu*m^2 * "]"),
     xlim = c(xlim_l, xlim_u),
     ylim = c(ylim_l, ylim_u),
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
points(d_input[-1], MSD_MFAIUQ[-1], pch = 1, col = "cyan", cex = 1.8, lwd = 1.8)
lines(d_input, truth, col = "black", lwd = 3)
polygon(c(d_input[-1], rev(d_input[-1])), c(MSD_lower[-1], rev(MSD_upper[-1])), col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)

# Microscopy image inset
img_raster = as.raster(img)
oldpar = par(no.readonly = TRUE)
par(fig = c(0.64, 0.92, 0.2, 0.59), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
rasterImage(img_raster, 0, 0, 1, 1)
scale_um = 50
field_um = pxsz * 512
bar_frac = scale_um/field_um
segments(0.62, 0.12, 0.62 + bar_frac, 0.12, col = "white", lwd = 5, lend = "butt")
par(oldpar)
legend("topleft", legend = legend_text,
       bty = "n", text.font = 1.8, cex = 2.2, inset = c(-0.1, -0.04))
