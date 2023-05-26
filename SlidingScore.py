# code from lab, not adjusted 

import numpy as np


nucs = {"A": 0, "C": 1, "G": 2, "T": 3} # this might be helpful

def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    for i in range(len(sequence)):
        score += pwm[nucs[sequence[i]], i]
    return score

def ReverseComplement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = []
    for i in sequence:
        if i == 'A':
            revcomp += 'T'
        elif i == 'T':
            revcomp += 'A'
        elif i == 'G':
            revcomp += 'C'
        else:
            revcomp += 'G'
    revcomp.reverse()
    revcomp = ("").join(revcomp)
    return revcomp


def FindMaxScore(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
    max_score = -1*np.inf
    # your code here
    n = pwm.shape[1]
    scores = [0]*(len(sequence)-n+1) # list of scores. scores[i] should give the score of the substring sequence[i:i+n]
    for i in range(len(scores)):
        reverse_comp = ReverseComplement(sequence[i:i+n])
        score_1 = ScoreSeq(pwm, sequence[i:i+n])
        score_2 = ScoreSeq(pwm, reverse_comp)
        if score_1 > score_2 and score_1 > max_score:
            max_score = score_1
        elif score_2 > score_1 and score_2 > max_score:
            max_score = score_2
    return max_score
