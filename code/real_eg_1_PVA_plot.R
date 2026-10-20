# Load required packages
library(bioimagetools)

# Plot the results using the provided records.
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

# Load the provided MF-AIUQ records
MFAIUQ_records = read.csv(paste0("../results/MF_AIUQ_records/MSD_", file_name, ".csv"), header = T)
d_input = as.matrix(MFAIUQ_records$d_input)[-1]
truth = as.matrix(MFAIUQ_records$truth)[-1]
MSD_MFAIUQ = as.matrix(MFAIUQ_records$n_par_AIUQ)[-1]
MSD_lower = as.matrix(MFAIUQ_records$lower)[-1]
MSD_upper = as.matrix(MFAIUQ_records$upper)[-1]

# Load the MPT results
tau_MPT = read.csv(paste0("../results/MPT_records/", file_name, "_tau.csv"), header = F)
tau_MPT = as.matrix(tau_MPT)
MSD_MPT = read.csv(paste0("../results/MPT_records/", file_name, "_MSD.csv"), header = F)
MSD_MPT = as.matrix(MSD_MPT)

# Load the DDM-UQ results
if (file_name == "4pct_PVA_100nm_25C_BF") {
  DDMUQ_records = read.csv(paste0("../results/DDM_UQ_records/DDMUQ_", file_name, ".csv"), header = T)
  tau_DDMUQ = as.matrix(DDMUQ_records[,1])
  MSD_DDMUQ = as.matrix(DDMUQ_records[,2])
} else {
  DDMUQ_records = read.csv(paste0("../results/DDM_UQ_records/DDMUQ_", file_name, ".csv"), header = F)
  tau_DDMUQ = as.matrix(DDMUQ_records[1,])
  MSD_DDMUQ = as.matrix(DDMUQ_records[2,])
}

# Load the inset intensity frame from the data/ folder
img = suppressWarnings(readTIF(paste0("../data/", file_name, "_inset.tif")))[,,1]

# Plot
if (scenario == "1 um") {
  ylim_l = min(MSD_lower) * 0.7
  legend_text = latex2exp::TeX("1 $\\mu$m")
} else if (scenario == "200 nm") {
  ylim_l = min(MSD_lower) * 0.1
  legend_text = "200 nm"
} else if (scenario == "100 nm") {
  ylim_l = min(MSD_lower) * 0.05
  legend_text = "100 nm"
} else if (scenario == "100 nm BF") {
  ylim_l = min(MSD_lower) * 0.05
  legend_text = "100 nm BF"
}

par(mar = c(4.1, 5.6, 1.5, 1.5), mgp = c(3.2, 1, 0), xpd = NA)
img = (img-min(img))/diff(range(img))
pxsz =0.2930 
xlim_l = 0.0309
xlim_u = max(d_input) + 2
ylim_u = max(range((c(truth, MSD_MFAIUQ, MSD_MPT, MSD_DDMUQ)))[2], max(MSD_upper)) * 1.1
plot(d_input, truth, log = "xy", type = "n",
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
points(d_input, MSD_MFAIUQ, pch = 1, col = "cyan", cex = 1.8, lwd = 1.8)
points(tau_DDMUQ, MSD_DDMUQ, pch = 2, col = "darkgreen", cex = 1.1, lwd = 1.1)
lines(tau_MPT, MSD_MPT, col = "#E76F51", lwd = 3)
lines(d_input, truth, col = "black", lwd = 3)
polygon(c(d_input, rev(d_input)), c(MSD_lower, rev(MSD_upper)), col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)
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



