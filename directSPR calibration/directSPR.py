import pydicom
import glob
import numpy as np
from math import sqrt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from collections import Counter

class LUT:
	def __init__(self, list_from, list_to):
		self.list_from = list_from
		self.list_to = list_to
		assert np.shape(self.list_from) == np.shape(self.list_to)
		self.interp = interp1d(self.list_from, self.list_to)

	def __call__(self, eval_at):
		return self.interp(eval_at)


class TissueType:
	def __init__(self, name, HUfrom, HUto):
		self.name = name
		self.HUfrom = HUfrom
		self.HUto = HUto
		self.SPRhistogram = { x : list() for x in range(HUfrom, HUto + 1) }

	def median(self, HU): # Median value of SPR with HU
		if self.getNumber(HU):
			return np.median(self.SPRhistogram[HU])
		else:
			return 0

	def getNumber(self, HU):
		return len(self.SPRhistogram[HU])

	def weight(self, HU): # All voxels with HU
		n = len(self.SPRhistogram[HU])
		if n>0:
			w = 1/sqrt(n)
			return w
		else:
			return 0

	def weights(self):
		yerr = list()
		for HU in sorted(self.SPRhistogram.keys()):
			yerr.append(self.weight(HU))
		return yerr


def map_rcs(mono, rho):
	if rho < 0.3:
		return 1
	else:
		mono_conv = (mono + 1000) / 1000.
		return mono_conv / rho

def map_spr(rsn, rho):
	return rsn * rho

######################################
# The following procedure is copied from section 2.2.2 in 10.1088/1361-6560/aaa1c9
# 1) rcs = Mono / Rho
# 2) Rho < 0.3 -> rcs = 1
# 3) rsn = rcs_to_rsn(rcs)
# 4) spr = rsn * Rho
######################################

######################################
# Table S3 in 10.1088/1361-6560/aaa1c9
rcs1 = [0, 0.9366, 0.9947, 1.0276, 1.1871, 1.6490, 2.2035, 9999]
rsn1_100 = [1.0286, 1.0286, 1.0028, 0.9955, 0.9949, 0.9495, 0.9070, 0.9070]
rsn1_200 = [1.0265, 1.0265, 1.0026, 0.9958, 0.9953, 0.9531, 0.9137, 0.9137]
rcs_to_rsn = LUT(rcs1, rsn1_200)
######################################

images_rho = dict()
images_mono = dict()
images_spr = dict()
images_rcs = dict()
images_rsn = dict()

# Dette er de gamle bildene
path_rho = "Sorted 60keV_RhoZ/DE Caput  2.0  Q30s  2 RhoZ Rho#1/*"
path_mono = "Sorted 60keV_RhoZ/DE Caput  2.0  Q30s  2 Monoenergetic Plus 60 keV#0/*"

fullpath_rho = glob.glob(path_rho)
fullpath_mono = glob.glob(path_mono)

print("Loading RHO dataset", end="")
for fullpath in fullpath_rho:
	print(".", end="")
	ds = pydicom.dcmread(fullpath)
	img = ds.pixel_array * ds.RescaleSlope + ds.RescaleIntercept
	n = ds.SliceLocation
	images_rho[n] = (img + 1000) / 1000
print("")

print("Loading MONO dataset", end="")
for fullpath in fullpath_mono:
	print(".", end="")
	ds = pydicom.dcmread(fullpath)
	img = ds.pixel_array * ds.RescaleSlope + ds.RescaleIntercept
	n = ds.SliceLocation
	images_mono[n] = img
print("")

tissues = [TissueType("Low Density", -900, -150), TissueType("Adipose", -125, -70), TissueType("MuscleBrain", -20, 70)]

print("Calculating SPR images")
for sl in images_rho.keys():
	print(f"Image at Slice Location {sl}: ")
	image_rho = images_rho[sl]
	image_mono = images_mono[sl]
	shape = np.shape(image_rho)
	rs_image_rho = image_rho.flatten()
	rs_image_mono = image_mono.flatten()

	print("Making RCS image...", end=" ")
	rs_image_rcs = np.array([map_rcs(i, j) for i,j in zip(rs_image_mono, rs_image_rho) ])
	images_rcs[sl] = np.reshape(rs_image_rcs, shape)
	print("Done!")

	print("Making RSN image...", end=" ")
	rs_image_rsn = [ rcs_to_rsn(i) for i in rs_image_rcs ]
	images_rsn[sl] = np.reshape(rs_image_rsn, shape)
	print("Done!")

	print("Making SPR image...", end=" ")
	rs_image_spr = np.array([map_spr(i, j) for i,j in zip(rs_image_rsn, rs_image_rho) ])
	images_spr[sl] = np.reshape(rs_image_spr, shape)

	for spr,hu in zip(rs_image_spr,rs_image_mono):
		for tissue in tissues:
			if tissue.HUfrom < hu < tissue.HUto:
				tissue.SPRhistogram[int(hu)].append(spr)

	print("Done!\n")
	break # Ta bort for å kjøre på hele datasettet, tar litt mer tid

fig, axs = plt.subplots(1, 5, figsize=(20,6), sharey=True)
plt.subplots_adjust(left = 0.04, right = 0.99, bottom=0, top=1, wspace=0.05)

k = list(images_mono.keys())[0]

axs[0].imshow(images_mono[k], cmap="gray")
axs[0].set_title("MonoCT 70 keV")

axs[1].imshow(images_rho[k], cmap="gray")
axs[1].set_title("Relative Electron Density")

axs[2].imshow(images_rcs[k], cmap="gray")
axs[2].set_title("Relative Cross Section")

axs[3].imshow(images_rsn[k], cmap="gray")
axs[3].set_title("Relative Stopping Number")

axs[4].imshow(images_spr[k], cmap="gray")
axs[4].set_title("Relative Stopping Power (directSPR)")

fig2, axs2 = plt.subplots(1, len(tissues), figsize=(20,6))
for idx, tissue in enumerate(tissues):
	x = list(tissue.SPRhistogram.keys())
	y = [ tissue.median(k) for k in x]
	yerr = [ tissue.weight(k) for k in x]

	while 0 in yerr:
		idx_0 = yerr.index(0)
		yerr.pop(idx_0)
		x.pop(idx_0)
		y.pop(idx_0)

	axs2[idx].errorbar(x, y, yerr=yerr)
	axs2[idx].set_title(tissue.name)

plt.show()
