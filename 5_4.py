# mass of attachments

# setup:
import numpy as np
import math

# constants:
massFuelTank = 817.75 #kg #placeholder
massRTG = 50 #kg
forceX = 1.05 * 9.80665 * massRTG #lateral static + dynamics loads ( 4.55 + 1) * g0 * mass of fuel tank
forceY = 1.05 * 9.80665 * massRTG #lateral static + dynamics loads ( 4.55 + 1) * g0 * mass of fuel tank
forceZ = 5.55 * 9.80665 * massRTG #longitudinal static + dynamics loads ( 4.55 + 1) * g0 * mass of fuelTank
#massLugBase = 10 #0.008966 #kg #mass of the lug calculated in wp4
totalForce = np.sqrt (forceX ** 2 + forceY ** 2 + forceZ ** 2) #total force for lug Base used in wp4
totalCompForce = 1000 #placeholder #total compressive force from the fuel tank

# variables
results = np.array(["Number of attachments", "Force ratio", "Attachment mass", "Total mass"]) #storing final results of iterations

# functions:

def massCalc (D_1,D_2,t_1,t_2, densityLug = 2810):
    attachmentVolume = (7.5 * D_2 * 2 * t_1) * t_2 * 6 * D_2 - D_2 ** 2 * np.pi * t_2 + 2 * (0.02 * t_1 * 6 * D_2 + (3 * D_2) ** 2 * t_1 * np.pi / 2 - t_1 / 4 * D_1 ** 2)
    attachmentMass = attachmentVolume * densityLug
    return attachmentMass

def massOfAttachments(totalForce, totalCompForce, numberOfAttachments):
    #calculating ratio of force per lug in this design compared to wp4
    forceRatio = 1/(numberOfAttachments * totalForce / totalCompForce)

    D_1 = 0.00265 * forceRatio
    D_2 = 0.00232 * forceRatio
    t_1 = 0.00418 * forceRatio
    t_2 = 0.00232 * forceRatio

    #calculating Mass of attachments as a ratio of forces in wp4 * size of lug in wp4
    attachmentMass = massCalc (D_1, D_2, t_1, t_2)

    #Total mass of attachment calculation
    totalAttachmentMass = numberOfAttachments * attachmentMass

    iteration = np.array([numberOfAttachments, forceRatio, attachmentMass, totalAttachmentMass])
    return iteration

# program
for i in range(1, 4):
    results = np.vstack((results, massOfAttachments(totalForce, totalCompForce, i)))

print (results)





