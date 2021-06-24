import glob
import numpy as np
import pandas as pd
from math import *
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pwlf


def sn(I, beta):
    return -(log(2 * 0.511 * beta**2 / (1 - beta**2)) - beta**2 - log(I))


mu_text = "Linear Attenuation [g/cm2]"
energy_text = "Energy [keV]"
rho_text = "Mass Density [g/cm3]"

photonEnergy = 70
protonEnergy = 20
protonMass = 938
Iw = 75

beta = sqrt(1 - protonMass**2 / (protonMass + protonEnergy)**2)

tissues = pd.read_csv("Data/Tissues/PMBaaa1c9_supplementarymaterial.csv", delimiter=",")
tissues.set_index("Insert", drop=True, inplace=True)
tissues["File"] = ["tissues_mohler2018"] * len(tissues)
tissues["RelativeElectronDensity"] = tissues["RED"]
tissues[rho_text] = tissues["RMD"]

tissues_bones = pd.read_csv("Data/Tissues/bones_white1987.csv", delimiter=",", header=0)
tissues_bones.drop(["Unnamed: 16"], axis=1, inplace=True)
tissues_bones["RelativeElectronDensity"] = tissues_bones["e/m3 x 10^26"] / 3346
tissues_bones.set_index("Tissue", drop=True, inplace=True)
tissues_bones["File"] = ["bones_white1987"] * len(tissues_bones)
tissues_bones[rho_text] = tissues_bones["Mass density (kg/m3)"] / 1000

tissues_tissues = pd.read_csv("Data/Tissues/tissues_woodard1986.csv", delimiter=",", header=0, engine="python")
tissues_tissues.drop(["Unnamed: 17", "Unnamed: 18"], axis=1, inplace=True)
tissues_tissues["RelativeElectronDensity"] = tissues_tissues["e/m3 x 10^26"] / 3346
tissues_tissues.set_index("Tissues", drop=True, inplace=True)
tissues_tissues.fillna(0, inplace=True)
tissues_tissues["File"] = ["tissues_woodard1986"] * len(tissues_tissues)
tissues_tissues[rho_text] = tissues_tissues["Mass density (kg/m3)"] / 1000

tissues = pd.concat([tissues, tissues_bones, tissues_tissues], join="inner")

tissues["lnI"] = np.zeros(len(tissues))
tissues["RelativeStoppingNumber"] = np.zeros(len(tissues))
tissues["RelativeCrossSection"] = np.zeros(len(tissues))
tissues["ElectronDensity"] = np.zeros(len(tissues))
tissues["MassAttenuation"] = np.zeros(len(tissues))

extraElementalData = pd.read_csv("Data/elementData.csv", delimiter=",")

elements = dict()
elementFiles = glob.glob("Data/Elements/TotalCrossSection/*.csv")
for elementFile in elementFiles:
    el = elementFile[:-4].split("\\")[1]
    elements[el] = pd.read_csv(elementFile, delimiter=" ", names=[energy_text, mu_text, "dummy"], usecols=[energy_text, mu_text])
    elements[el][energy_text] *= 1000

    # Interpolate the cross section if the input photon energy is not in the list
    if photonEnergy not in list(elements[el][energy_text]):
        elements[el] = elements[el].append(pd.Series({energy_text: photonEnergy, mu_text: np.nan}), ignore_index=True)
        elements[el].set_index(energy_text, drop=True, inplace=True)
        elements[el].sort_values(energy_text, inplace=True)
        elements[el].interpolate(inplace=True, method="slinear")
    else:
        elements[el].set_index(energy_text, drop=True, inplace=True)


# Calculate cross section for water
MassAttenuationWater = 0
for element in tissues.loc["Water"].index:
    if element not in elements:
        continue

    wi = float(tissues.at["Water", element]) / 100
    MassAttenuationForElement = elements[element].loc[float(photonEnergy)][mu_text]
    MassAttenuationWater += MassAttenuationForElement * wi

LinearAttenuationWater = MassAttenuationWater  # since rho = 1
CrossSectionWater = LinearAttenuationWater  # since n = 1

# Electron density normalization
for tissue in tissues.index:
    for element in tissues.loc[tissue].index:
        if element not in elements:
            continue

        wi = float(tissues.at[tissue, element]) / 100
        Zi = float(extraElementalData[extraElementalData["Element"] == element]["Z"])
        Ai = float(extraElementalData[extraElementalData["Element"] == element]["A"])
        tissues.at[tissue, "ElectronDensity"] += wi * (Zi / Ai)

# Calculate relative cross section
for tissue in tissues.index:
    RelativeElectronDensity = float(tissues.at[tissue, "RelativeElectronDensity"])
    MassDensity = float(tissues.at[tissue, rho_text])
    for element in tissues.loc[tissue].index:
        if element not in elements:
            continue

        wi = float(tissues.at[tissue, element]) / 100
        MassAttenuationForElement = elements[element].loc[float(photonEnergy)][mu_text]
        tissues.at[tissue, "MassAttenuation"] += MassAttenuationForElement * wi

    LinearAttenuation = tissues.at[tissue, "MassAttenuation"] * MassDensity
    RelativeLinearAttenuation = LinearAttenuation / LinearAttenuationWater
    RelativeCrossSection = RelativeLinearAttenuation / RelativeElectronDensity
    tissues.at[tissue, "RelativeCrossSection"] = RelativeCrossSection

# Calculate relative stopping number
for tissue in tissues.index:
    MassDensity = float(tissues.at[tissue, rho_text])
    RelativeElectronDensity = float(tissues.at[tissue, "RelativeElectronDensity"])
    for element in tissues.loc[tissue].index:
        if element not in elements:
            continue

        wi = float(tissues.at[tissue, element]) / 100
        Zi = float(extraElementalData[extraElementalData["Element"] == element]["Z"])
        Ai = float(extraElementalData[extraElementalData["Element"] == element]["A"])
        Ii = float(extraElementalData[extraElementalData["Element"] == element]["I"])

        RelativeStoppingNumberForElement = sn(Ii, beta) / sn(Iw, beta)

        # nu_i in Eq (3) of PMB 61 N286
        nu = (wi * (Zi / Ai)) / float(tissues.at[tissue, "ElectronDensity"])

        tissues.at[tissue, "RelativeStoppingNumber"] += RelativeStoppingNumberForElement * nu

# Print & plot
for tissue in tissues.index:
    RelativeElectronDensity = float(tissues.at[tissue, "RelativeElectronDensity"])
    MassDensity = float(tissues.at[tissue, rho_text])
    rcs = tissues.at[tissue, "RelativeCrossSection"]
    rsn = tissues.at[tissue, "RelativeStoppingNumber"]
    rsp = tissues.at[tissue, "RelativeStoppingNumber"] * RelativeElectronDensity
    kind = tissues.at[tissue, "File"]

    print(f"{kind:<20} / {tissue:<45}\t: RED = {RelativeElectronDensity:.3f}; RCS {photonEnergy} keV = {rcs:.3f}; "
          f"RSN {protonEnergy} MeV = {rsn:.3f}; RSP = {rsp:.3f}")

mkr_color = {'tissues_mohler2018': 'red', 'bones_white1987': 'blue', 'tissues_woodard1986': 'green'}

rsp_hunemohr = {'Lung': 0.444, 'Adipose': 0.943, 'Breast': 0.983, 'Water': 1.0, 'Solid water': 1.001, 'Muscle': 1.033,
                'Brain': 1.064, 'Liver': 1.073, 'Inner bone': 1.099, 'B-200': 1.108, 'CB2-30%': 1.263, 'CB2-50%': 1.426, 'Cortical bone': 1.612}

for kind in mkr_color:
    d = tissues[tissues.File == kind]
    plt.scatter(d.RelativeCrossSection, d.RelativeStoppingNumber, marker="o", s=40, alpha=0.4, c=mkr_color[kind], label=kind)
    plt.xlabel("Relative Cross Section")
    plt.ylabel("Relative Stopping Number")
    plt.legend()

x0 = [0, 1, 1.2, tissues.RelativeCrossSection.max() * 2]
my_pwlf = pwlf.PiecewiseLinFit(tissues.RelativeCrossSection, tissues.RelativeStoppingNumber)
my_pwlf.fit_with_breaks(x0)

xHat = np.linspace(x0[0], x0[-1], 5000)
yHat = my_pwlf.predict(xHat)
print(my_pwlf.calc_slopes())
print(my_pwlf.intercepts)

plt.plot(xHat, yHat, "-")

plt.xlim([0.85, 1.75])
plt.ylim([0.925, 1.15])

plt.show()

print("\nRelative Cross Section -> Relative Stopping Number LUT.\nBreakpoints:")
for x in x0:
    print(f"RCS = {x:.4f} -> RSN = {my_pwlf.predict(x)[0]:.4f}")

"""
# Check RSP accuracy where data is available (Hünemohr for standard tissues)
plt.figure()

d = tissues[tissues.File == "tissues_mohler2018"]
x = np.array(list(d.RelativeStoppingNumber * d.RelativeElectronDensity))
y = np.array([rsp_hunemohr[k] for k in tissues.index if k in rsp_hunemohr])
plt.scatter(x, (x - y) / y * 100)
plt.xlabel("My Relative Stopping Power")
plt.ylabel("RSP error compared to Hünemohr [%]")

plt.show()
"""
