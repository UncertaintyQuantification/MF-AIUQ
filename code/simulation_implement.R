# Load required package
library(AIUQ)

# Select one scenario:
# "Slow BM", "Fast BM", "Subdiffusive FBM",
# "Superdiffusive FBM", "OU", or "OU-FBM"
scenario = "Slow BM"

# Simulation and specific settings for plots
if (scenario == "Slow BM") {
  simul_model_name = "BM"
  sigma_bm = 0.1
  legend_text = "Slow BM"
  ylim_l = 0.01
  } else if (scenario == "Fast BM") {
  simul_model_name = "BM"
  sigma_bm = 1
  legend_text = "Fast BM"
  ylim_l = 1
} else if (scenario == "Subdiffusive FBM") {
  simul_model_name = "FBM"
  sigma_fbm = 2
  H = 0.3
  legend_text = "Subdiffusive FBM"
  ylim_l = 1
} else if (scenario == "Superdiffusive FBM") {
  simul_model_name = "FBM"
  sigma_fbm = 0.5
  H = 0.6
  legend_text = "Superdiffusive FBM"
  ylim_l = 0.1
} else if (scenario == "OU") {
  simul_model_name = "OU"
  sigma_ou = 4
  rho = 0.95
  legend_text = "OU"
  ylim_l = 1
} else if (scenario == "OU-FBM") {
  simul_model_name = "OU+FBM"
  sigma_ou = 1.5
  rho = 0.85
  sigma_fbm = 1
  H = 0.225
  legend_text = "OU-FBM"
  ylim_l = 1
}


set.seed(1)
M = 100
frame_size = 500
len = 500
I0 = 20
noise_type = "gaussian"
pos0 = matrix(0.25*frame_size + 0.5*frame_size*runif(2*M), nrow = M, ncol = 2)
# Simulate the selected scenario
if (scenario == "Slow BM" || scenario == "Fast BM") {
  simulation_list = AIUQ::simulation(sz = frame_size, len_t = len, model_name = simul_model_name, noise = noise_type, I0 = I0, pos0 = pos0, M = M, sigma_bm = sigma_bm)
} else if (scenario == "Subdiffusive FBM" || scenario == "Superdiffusive FBM") {
  simulation_list = AIUQ::simulation(sz = frame_size, len_t = len, model_name = simul_model_name, noise = noise_type, I0 = I0, pos0 = pos0, M = M, sigma_fbm = sigma_fbm, H = H)
} else if (scenario == "OU") {
  simulation_list = AIUQ::simulation(sz = frame_size, len_t = len, model_name = simul_model_name, noise = noise_type, I0 = I0, pos0 = pos0, M = M, sigma_ou = sigma_ou, rho = rho)
} else if (scenario == "OU-FBM") {
  simulation_list = AIUQ::simulation(sz = frame_size, len_t = len, model_name = simul_model_name, noise = noise_type, I0 = I0, pos0 = pos0, M = M, sigma_ou = sigma_ou, rho = rho, sigma_fbm = sigma_fbm, H = H)
}

# Estimate the MSD using MF-AIUQ
# Set uncertainty = FALSE if confidence intervals are not required
sam_mf = SAM(sim_object = simulation_list, model_name = simulation_list@model_name,
             model_free_estimation = TRUE, uncertainty = TRUE)
# MF-AIUQ results
d_input = sam_mf@d_input
MSD_MFAIUQ = sam_mf@msd_est
MSD_lower = sam_mf@msd_lower
MSD_upper = sam_mf@msd_upper
truth = simulation_list@theor_msd

# Plot
par(mar = c(3.7, 5.6, 1.5, 1.5), mgp = c(2.8, 1, 0), xpd = NA)
xlim_l = 1
xlim_u = max(d_input[-1]) + 60
ylim_l = ylim_l
ylim_u = max(range((c(truth[-1], MSD_MFAIUQ[-1])))[2], max(MSD_upper[-1])) * 1.1
plot(d_input[-1], truth[-1], log = "xy", type = "n",
     xlab = expression(Delta*t), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">"),
     xlim = c(xlim_l, xlim_u),
     ylim = c(ylim_l, ylim_u),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", cex.lab = 2.5)
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
lines(d_input[-1], truth[-1], col = "black", lwd = 3)
polygon(c(d_input[-1], rev(d_input[-1])),
        c(MSD_lower[-1], rev(MSD_upper[-1])),
        col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)
legend("bottomright", legend = legend_text,
       bty = "n", text.font = 1.8, cex = 2.2, inset = c(-0.03, -0.04))
