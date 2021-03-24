
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
 
    
def read_substitution_matrix(file_substitution_matrix):
    """
    Implement reading the scores file.
    It can be stored as a dictionary of example:
    scores[("A", "R")] = -1
    """
    scores = {}
    with open(file_substitution_matrix) as f:
        #to get only the lines belonging to scoring matrix
        lines = [line.split() for line in f if not line.startswith('#') and not line.startswith('\n')]
    
    #coverting list into a dictionary
    for i in range(21):
        for j in range(1, 22):
            scores[(lines[j][0], lines[0][i])] = int(lines[j][i+1])
    
    return scores
    
    
def init_matrix_d(seq_1, seq_2, cost_gap_open, cost_gap_extend):
    """
    Implement initialization of the matrix D
    """
    matrix_d = []
    
    while len(matrix_d) < len(seq_1)+1:
        matrix_d.append([])
        while len(matrix_d[-1]) < len(seq_2)+1:
            matrix_d[-1].append(0)
            
    for i in range(1, len(matrix_d[0])):
        matrix_d[0][i] = cost_gap_open + cost_gap_extend*(i)
    for j in range(1, len(matrix_d)):
        matrix_d[j][0] = cost_gap_open + cost_gap_extend*(j)
        
    return matrix_d
    

def init_matrix_p(seq_1, seq_2):
    """
    Implement initialization of the matrix P
    """
    matrix_p = []
    
    while len(matrix_p) < len(seq_1)+1:
        matrix_p.append([])
        while len(matrix_p[-1]) < len(seq_2)+1:
            matrix_p[-1].append(0)
    
    for i in range(1, len(matrix_p[0])):
        matrix_p[0][i] = -float('Inf')
    for j in range(1, len(matrix_p)):
        matrix_p[j][0] = float('NaN')       
    
    return matrix_p
    
    
def init_matrix_q(seq_1, seq_2):
    """
    Implement initialization of the matrix Q
    """
    matrix_q = []
    
    while len(matrix_q) < len(seq_1)+1:
        matrix_q.append([])
        while len(matrix_q[-1]) < len(seq_2)+1:
            matrix_q[-1].append(0)
    
    for i in range(1, len(matrix_q[0])):
        matrix_q[0][i] = float('NaN')
    for j in range(1, len(matrix_q)):
        matrix_q[j][0] = -float('Inf')
        
    return matrix_q
