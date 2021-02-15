
def score_of_alignment(align_seq1, align_seq2, cost_gap_open, 
                       cost_gap_extension, substitutions=None):
    """
    A nice helper function which computes the score of the given alignment.
    This is only used for the self check.
    Input example:
    --CG
    AACA
    """
    score = 0
    i = 0
    
    if substitutions == None:
        while i < len(align_seq1):
            if align_seq1[i] == align_seq2[i]:
                score += 1
                i += 1
            elif align_seq1[i] == '-':
                extension = 1
                for j in range(i+1, len(align_seq1)):
                    if align_seq1[j] == '-':
                        extension += 1
                        i += 1
                    else:
                        break
                score += cost_gap_open + cost_gap_extension * extension
                i += 1
            elif align_seq2[i] == '-':
                extension = 1
                for j in range(i+1, len(align_seq2)):
                    if align_seq2[j] == '-':
                        extension += 1
                        i += 1
                    else:
                        break
                score += cost_gap_open + cost_gap_extension * extension
                i += 1
            else:
                score = score - 1
                i += 1
                
    else:
        while i < len(align_seq1):
            if align_seq1[i] == '-':
                extension = 1
                for j in range(i+1, len(align_seq1)):
                    if align_seq1[j] == '-':
                        extension += 1
                        i += 1
                    else:
                        break
                score += cost_gap_open + cost_gap_extension * extension
                i += 1
            elif align_seq2[i] == '-':
                extension = 1
                for j in range(i+1, len(align_seq2)):
                    if align_seq2[j] == '-':
                        extension += 1
                        i += 1
                    else:
                        break
                score += cost_gap_open + cost_gap_extension * extension
                i += 1
            else:
                score += substitutions[(align_seq1[i], align_seq2[i])]
                i += 1               
    
    return score
    

def read_fasta_file(fasta_file):
    "Implement reading the fasta file"
    with open(fasta_file) as f:
        lines = [line.rstrip() for line in f]
        #remove the first line which is not part of the sequence
        sequence = lines[1]
    return sequence
