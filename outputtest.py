import seqlogo
import pickle
import numpy as np
import logomaker as lm
import pandas as pd
import matplotlib.pyplot as plt

def output(pwm):
    # make seqlogo PWM object
    #seq_pwm = seqlogo.Pwm(pwm)
    # Convert to ppm needed for plotting
    #seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
    #seqlogo.seqlogo(seq_ppm, ic_scale = True, format = 'png', size = 'medium')

    # different package: LogoMaker
    rows = []
    for i in range(0, len(pwm)):
        rows.append(str(i + 1))
    pwm = pd.DataFrame(data=pwm, columns=["A","C","G","T"])
    logo = lm.Logo(pwm, font_name='Arial')
    return logo
    
motifs = pickle.load(open("Motif_Find/Jaspar.p", "rb"))
pwm = np.array(motifs["MA0139.1"][0]).transpose()
logo = output(pwm)
plt.savefig("logo", format='png')
# output(pwm)

# # Sox2
# pfm_1 = np.array([
# [35.0,92.0,146.0,226.0],
# [28.0,281.0,71.0,119.0],
# [25.0,272.0,8.0,194.0],
# [219.0,6.0,4.0,270.0],
# [1.0,11.0,5.0,482.0],
# [7.0,1.0,2.0,489.0],
# [39.0,76.0,372.0,12.0],
# [52.0,1.0,12.0,434.0],
# [113.0,124.0,78.0,184.0],
# [340.0,39.0,33.0,87.0],
# [36.0,5.0,19.0,439.0],
# [19.0,29.0,375.0,76.0],
# [24.0,317.0,59.0,99.0],
# [279.0,51.0,17.0,152.0],
# [277.0,47.0,72.0,103.0],
# [386.0,19.0,38.0,56.0]
# ]).transpose()

# # Gata4
# pfm_2 = np.array([
# [184.0,64.0,156.0,96.0],
# [94.0,215.0,142.0,49.0],
# [332.0,4.0,8.0,156.0],
# [2.0,0.0,498.0,0.0],
# [490.0,2.0,4.0,4.0],
# [4.0,11.0,20.0,465.0],
# [481.0,2.0,2.0,15.0],
# [459.0,1.0,35.0,5.0],
# [47.0,91.0,349.0,13.0],
# [285.0,71.0,128.0,16.0],
# ]).transpose()

# # HNF4A
# pfm_3 = np.array([
#     [141.0,127.0,177.0,55.0],
# [228.0,16.0,199.0,57.0],
# [36.0,9.0,424.0,31.0],
# [91.0,35.0,163.0,211.0],
# [75.0,154.0,139.0,132.0],
# [4.0,463.0,5.0,28.0],
# [475.0,11.0,8.0,6.0],
# [423.0,4.0,72.0,1.0],
# [461.0,0.0,34.0,5.0],
# [4.0,3.0,484.0,9.0],
# [8.0,4.0,191.0,297.0],
# [21.0,248.0,64.0,167.0],
# [9.0,418.0,8.0,65.0],
# [351.0,32.0,65.0,52.0]
# ]).transpose()

# # Compute the pwms, which we'll use below
# PWMList = []
# pwm_names = ["SOX2","GATA4","HNF4A"]

# # Plot the sequence logos
# for pfm in [pfm_1, pfm_2, pfm_3]:
#     seq_pfm = seqlogo.Pfm(pfm/np.sum(pfm, 0)[0]) # normalize to probabilities rather than counts
#     seq_ppm = seqlogo.Ppm(seqlogo.pfm2ppm(seq_pfm))
#     PWMList.append(np.array(seqlogo.ppm2pwm(seq_ppm)).transpose())
#     seqlogo.seqlogo(seq_ppm, ic_scale = True, format = 'png', size = 'medium')

# print("Example PWM")
# PWMList[0]