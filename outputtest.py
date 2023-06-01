import seqlogo
import pickle

def output(pwm):
    # make seqlogo PWM object
    seq_pwm = seqlogo.Pwm(pwm)
    # Convert to ppm needed for plotting
    seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
    seqlogo.seqlogo(seq_ppm, ic_scale = True, format = 'png', size = 'medium')
    
motifs = pickle.load(open("Motif_Find/Jaspar.p", "rb"))
pwm = motifs["MA0139.1"][0]
print(pwm)
output(pwm)