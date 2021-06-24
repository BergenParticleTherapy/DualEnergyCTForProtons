import glob
import numpy as np
import pandas as pd
import pydicom
from math import *
from matplotlib import pyplot as plt


def ROI(ar, circle, r):
    ROI = np.zeros(np.shape(ar), dtype=np.bool_)
    x0, y0 = circle

    for y in range(floor(y0 - r), ceil(y0 + r) + 1):
        for x in range(floor(x0 - r), ceil(x0 + r) + 1):
            if (x - x0)**2 + (y - y0)**2 < r**2:
                ROI[y, x] = True
    return ROI


def findEffectiveEnergy(tissues, k, HU, std):
    HU_error = dict()
    HU_error_upper = dict()
    HU_error_lower = dict()

    for energy in energies:
        calcHU = float(tissues.at[k, f"HounsfieldUnit_{energy}"])
        HU_error[energy] = abs(HU - calcHU)
        HU_error_upper[energy] = abs(HU + std / 2 - calcHU)
        HU_error_lower[energy] = abs(HU - std / 2 - calcHU)

    bestfit = min(HU_error, key=HU_error.get)
    bestfit_upper = min(HU_error_upper, key=HU_error_upper.get)
    bestfit_lower = min(HU_error_lower, key=HU_error_lower.get)

    error = abs(bestfit_upper - bestfit_lower)
    if not error:
        error = 1

    return bestfit, error


mu_text = "Linear Attenuation [g/cm2]"
energy_text = "Energy [keV]"
rho_text = "Mass Density [g/cm3]"

images_80 = "Data/Bildeserier/import/DE_Abdomen  5.0  D40f  A_80kV/*"
images_140 = "Data/Bildeserier/import/DE_Abdomen  5.0  D40f  B_Sn140kV/*"
imageNo = 7
energies = range(10, 500)

ds_80 = pydicom.dcmread(glob.glob(images_80)[imageNo])
ds_140 = pydicom.dcmread(glob.glob(images_140)[imageNo])

tissues = pd.read_csv("Data/Tissues/PMBaaa1c9_supplementarymaterial.csv", delimiter=",")
tissues.set_index("Insert", drop=True, inplace=True)
tissues[rho_text] = tissues["RMD"]

for energy in energies:
    tissues[f"MassAttenuation_{energy}"] = np.zeros(len(tissues))
    tissues[f"HounsfieldUnit_{energy}"] = np.zeros(len(tissues))

extraElementalData = pd.read_csv("Data/elementData.csv", delimiter=",")

elements = dict()
elementFiles = glob.glob("Data/Elements/TotalCrossSection/*.csv")
for elementFile in elementFiles:
    el = elementFile[:-4].split("\\")[1]
    elements[el] = pd.read_csv(elementFile, delimiter=" ", names=[energy_text, mu_text, "dummy"], usecols=[energy_text, mu_text])
    elements[el][energy_text] *= 1000

    # Add energies
    for energy in energies:
        if energy not in list(elements[el][energy_text]):
            elements[el] = elements[el].append(pd.Series({energy_text: energy, mu_text: np.nan}), ignore_index=True)

    elements[el].set_index(energy_text, drop=True, inplace=True)
    elements[el].sort_values(energy_text, inplace=True)

    # Logarithmic interpolation for higher accuracy (cf. XCOM data visually)
    elements[el][mu_text] = np.log10(elements[el][mu_text])
    elements[el].loc[energy] = np.log10(elements[el].loc[energy])
    elements[el].interpolate(inplace=True, method="slinear")
    elements[el][mu_text] = np.power(10, elements[el][mu_text])
    elements[el].loc[energy] = np.power(10, elements[el].loc[energy])


# Calculate mass attenuation section for water
for element in tissues.loc["Water"].index:
    if element not in elements:
        continue

    for energy in energies:
        wi = float(tissues.at["Water", element]) / 100
        MassAttenuationForElement = elements[element].loc[float(energy)][mu_text]
        tissues.at["Water", f"MassAttenuation_{energy}"] += MassAttenuationForElement * wi


# Calculate HU @ different energies
for tissue in tissues.index:
    if tissue == "Water":
        continue

    MassDensity = float(tissues.at[tissue, rho_text])
    for element in tissues.loc[tissue].index:
        if element not in elements:
            continue

        for energy in energies:
            wi = float(tissues.at[tissue, element]) / 100
            MassAttenuationForElement = elements[element].loc[float(energy)][mu_text]
            tissues.at[tissue, f"MassAttenuation_{energy}"] += MassAttenuationForElement * wi

    for energy in energies:
        LinearAttenuation = tissues.at[tissue, f"MassAttenuation_{energy}"] * MassDensity
        RelativeLinearAttenuation = LinearAttenuation / tissues.at["Water", f"MassAttenuation_{energy}"]
        tissues.at[tissue, f"HounsfieldUnit_{energy}"] = (RelativeLinearAttenuation - 1) * 1000


# Load gammex images and find mean HU values in ROIs
img_80 = ds_80.pixel_array * ds_80.RescaleSlope + ds_80.RescaleIntercept
img_140 = ds_140.pixel_array * ds_140.RescaleSlope + ds_140.RescaleIntercept

alpha = 1.496
img_mono = alpha * img_140 + (1 - alpha) * img_80

radius = 15
circles = {"CB2-30%": (189, 105), "Lung": (313, 105), "Breast": (191, 194),
           "SW1": (250, 169), "Adipose": (310, 195), "Inner bone": (399, 194),
           "Water": (166, 255), "SW2": (335, 255), "Brain": (101, 317),
           "CB2-50%": (190, 314), "Muscle": (251, 340), "B-200": (310, 315),
           "SW3": (398, 317), "Liver": (187, 405), "Cortical bone": (313, 405)}


plt.imshow(img_80, cmap="gray")
ax_img = plt.gca()

SW1 = np.mean(img_80[ROI(img_80, circles["SW1"], radius)])
SW2 = np.mean(img_80[ROI(img_80, circles["SW2"], radius)])
SW3 = np.mean(img_80[ROI(img_80, circles["SW3"], radius)])

SW1 = (SW1 + SW2) / 2

# Correct for beam hardening / ROI radial position in phantom
sw_correct = {'CB2-30%': SW3, 'Lung': SW3, 'Inner bone': SW3, 'Cortical bone': SW3, 'Liver': SW3, 'Brain': SW3,
              'Breast': SW1, 'Adipose': SW1, 'B-200': SW1, 'Muscle': SW1, 'CB2-50%': SW1, 'Breast': SW1, 'Water': SW1}

bestfits_80 = dict()
bestfits_140 = dict()
bestfits_mono = dict()

calcHUerror_80 = dict()
for energy in energies:
    calcHUerror_80[energy] = 0

calcHUerror_140 = dict()
for energy in energies:
    calcHUerror_140[energy] = 0

calcHUerror_mono = dict()
for energy in energies:
    calcHUerror_mono[energy] = 0

for tissue, v in circles.items():
    if tissue in ["SW1", "SW2", "SW3"]:
        continue

    c = plt.Circle(v, radius, alpha=0.5, color="red")
    ax_img.add_artist(c)

    filter = ROI(img_80, v, radius)
    ROIi = img_80[filter]
    mean = np.mean(ROIi)
    std = np.std(ROIi)
    mean_corrected = mean - sw_correct[tissue]
    effective_energy, effective_energy_error = findEffectiveEnergy(tissues, tissue, mean_corrected, std)

    # Best effective energy per insert
    bestfits_80[tissue] = (effective_energy, effective_energy_error)

    # Best single effective energy for all inserts
    for energy in energies:
        calcHUerror_80[energy] += (mean_corrected - float(tissues.at[tissue, f"HounsfieldUnit_{energy}"])) ** 2

    print(f"{tissue} (80 kVp) has an HU of {mean:.2f} +- {std:.2f} -> effective energy is {effective_energy} keV")


# New calculation: RMS of all modules per effective energy
bestFit = min(calcHUerror_80, key=calcHUerror_80.get)
X = list(calcHUerror_80.keys())
Y = [sqrt(v) for v in calcHUerror_80.values()]

plt.figure()
plt.plot(X, Y)
plt.title("RMS fit of effective energies (80 kVp)")
plt.xlabel("Effective energy [keV]")
plt.ylabel("RMS")
print(f"Overall best fit @ 80 kVp is {bestFit} keV.")
##########

fig, axs = plt.subplots(1, 2, figsize=(15, 8))

labels = bestfits_80.keys()
x_pos = np.arange(len(labels))
y = [k[0] for k in bestfits_80.values()]
yerr = [k[1] for k in bestfits_80.values()]
average = np.mean(y)

axs[0].bar(x_pos, y, yerr=yerr, align='center', alpha=0.5, ecolor='black', capsize=10)
axs[0].set_title(f"Effective energy per material at 80 kVp. Average = {average:.1f} keV.")
axs[0].set_xticks(x_pos)
plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=45)
axs[0].set_xticklabels(labels)
axs[0].plot([0, x_pos[-1]], [average] * 2, '-', color='red')

SW1 = np.mean(img_140[ROI(img_140, circles["SW1"], radius)])
SW2 = np.mean(img_140[ROI(img_140, circles["SW2"], radius)])
SW3 = np.mean(img_140[ROI(img_140, circles["SW3"], radius)])

SW1 = (SW1 + SW2) / 2

sw_correct = {'CB2-30%': SW3, 'Lung': SW3, 'Inner bone': SW3, 'Cortical bone': SW3, 'Liver': SW3, 'Brain': SW3,
              'Breast': SW1, 'Adipose': SW1, 'B-200': SW1, 'Muscle': SW1, 'CB2-50%': SW1, 'Breast': SW1, 'Water': SW1}

for tissue, v in circles.items():
    if tissue in ["SW1", "SW2", "SW3"]:
        continue

    filter = ROI(img_140, v, radius)
    ROIi = img_140[filter]
    mean = np.mean(ROIi)
    std = np.std(ROIi)
    mean_corrected = mean - sw_correct[tissue]
    effective_energy, effective_energy_error = findEffectiveEnergy(tissues, tissue, mean_corrected, std)
    bestfits_140[tissue] = (effective_energy, effective_energy_error)

    # Best single effective energy for all inserts
    for energy in energies:
        calcHUerror_140[energy] += (mean_corrected - float(tissues.at[tissue, f"HounsfieldUnit_{energy}"])) ** 2

    print(f"{tissue} (140 kVp) has an HU of {mean:.2f} +- {std:.2f} -> effective energy is {effective_energy} keV")

labels = bestfits_140.keys()
x_pos = np.arange(len(labels))
y = [k[0] for k in bestfits_140.values()]
yerr = [k[1] for k in bestfits_140.values()]
average = np.mean(y)  # Maybe weigh this somehow ? 1/yerr^2 doesn't seem right (np.average(y, weights = my_weights))

axs[1].bar(x_pos, y, yerr=yerr, align='center', alpha=0.5, ecolor='black', capsize=10)
axs[1].set_title(f"Effective energy per material at 140 kVp. Average = {average:.1f} keV.")
axs[1].set_xticks(x_pos)
plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=45)
axs[1].set_xticklabels(labels)
axs[1].plot([0, x_pos[-1]], [average] * 2, '-', color='red')

axs[0].set_ylim(0, 210)
axs[1].set_ylim(0, 210)

# New calculation: RMS of all modules per effective energy
bestFit = min(calcHUerror_140, key=calcHUerror_140.get)
X = list(calcHUerror_140.keys())
Y = [sqrt(v) for v in calcHUerror_140.values()]

plt.figure()
plt.plot(X, Y)
plt.title("RMS fit of effective energies (140 kVp)")
plt.xlabel("Effective energy [keV]")
plt.ylabel("RMS")
print(f"Overall best fit @ 140 kVp is {bestFit} keV.")
##########


# MONO

SW1 = np.mean(img_mono[ROI(img_mono, circles["SW1"], radius)])
SW2 = np.mean(img_mono[ROI(img_mono, circles["SW2"], radius)])
SW3 = np.mean(img_mono[ROI(img_mono, circles["SW3"], radius)])

SW1 = (SW1 + SW2) / 2

sw_correct = {'CB2-30%': SW3, 'Lung': SW3, 'Inner bone': SW3, 'Cortical bone': SW3, 'Liver': SW3, 'Brain': SW3,
              'Breast': SW1, 'Adipose': SW1, 'B-200': SW1, 'Muscle': SW1, 'CB2-50%': SW1, 'Breast': SW1, 'Water': SW1}

for tissue, v in circles.items():
    if tissue in ["SW1", "SW2", "SW3"]:
        continue

    filter = ROI(img_mono, v, radius)
    ROIi = img_mono[filter]
    mean = np.mean(ROIi)
    std = np.std(ROIi)
    mean_corrected = mean - sw_correct[tissue]
    effective_energy, effective_energy_error = findEffectiveEnergy(tissues, tissue, mean_corrected, std)
    bestfits_mono[tissue] = (effective_energy, effective_energy_error)

    # Best single effective energy for all inserts
    for energy in energies:
        calcHUerror_mono[energy] += (mean_corrected - float(tissues.at[tissue, f"HounsfieldUnit_{energy}"])) ** 2

    print(f"{tissue} (mono) has an HU of {mean:.2f} +- {std:.2f} -> effective energy is {effective_energy} keV")

labels = bestfits_mono.keys()
x_pos = np.arange(len(labels))
y = [k[0] for k in bestfits_mono.values()]
yerr = [k[1] for k in bestfits_mono.values()]
average = np.mean(y)  # Maybe weigh this somehow ? 1/yerr^2 doesn't seem right (np.average(y, weights = my_weights))

# New calculation: RMS of all modules per effective energy
bestFit = min(calcHUerror_mono, key=calcHUerror_mono.get)
X = list(calcHUerror_mono.keys())
Y = [sqrt(v) for v in calcHUerror_mono.values()]

plt.figure()
plt.plot(X, Y)
plt.title("RMS fit of effective energies (mono kVp)")
plt.xlabel("Effective energy [keV]")
plt.ylabel("RMS")
print(f"Overall best fit @ mono kVp is {bestFit} keV.")
##########

plt.show()
