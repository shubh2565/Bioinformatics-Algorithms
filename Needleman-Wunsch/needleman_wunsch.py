import numpy as np

#main function
def needleman_wunsch(fasta_file_1, fasta_file_2, file_substitution_matrix, cost_gap_open):
    "Put your code here"
    
    seq1 = read_fasta_file(fasta_file_1)
    seq2 = read_fasta_file(fasta_file_2)
    scores = read_substitution_matrix(file_substitution_matrix)
    
    alignment_matrix = calculate_matrix(seq1, seq2, cost_gap_open, scores)
    trace = traceback(alignment_matrix, seq1, seq2, cost_gap_open, scores)
    
    aligned_seq_1, aligned_seq_2 = alingment_build(trace, seq1, seq2)
    alignment_score = alignment_matrix[-1][-1]
    
    return alignment_score, aligned_seq_1, aligned_seq_2

#function to read the sequence from a fasta file
def read_fasta_file(fasta_file):
    "Implement reading the fasta file"
    with open(fasta_file) as f:
        lines = [line.rstrip() for line in f]
        #remove the first line which is not part of the sequence
        sequence = lines[1]
    return sequence

#function to get the scores from pam or blosum file in a dictionary
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

#intializationo of matrix
def init_matrix(sequence1, sequence2, gap_cost):
    """
    Implement initialization of the matrix.
    Make sure you picked the right dimention and correctly initilized the first row and the first column.
    """
    matrix = np.zeros((len(sequence1)+1, len(sequence2)+1))
    matrix[:,0] = [i*gap_cost for i in range(np.shape(matrix)[0])]
    matrix[0,:] = [i*gap_cost for i in range(np.shape(matrix)[1])]
    
    return matrix

#compuation of new value
def new_value_computation(char_seq1, char_seq2, gap_cost, substitution_matrix, diag_val, top_val, left_val):
    """
    Implement the computation of the value in the new cell.
    In this function we assume that we want to compute the value in the new cell in the matrix.
    Assume that the values "to the left", "to the top" and "top left" are already computed and provided
    as the input to the function. Also we know what characters in both sequences correspond to the given cell.
    """
    diag = diag_val + substitution_matrix[(char_seq1, char_seq2)]
    top = top_val + gap_cost
    left = left_val + gap_cost
    cell_value = max(diag, top, left)
    
    return cell_value

def calculate_matrix(sequence1, sequence2, gap_opening_cost, substitution_cost):
    """
    Implement the step of complete computation of the matrix
    First initialize the matrix then fill it in from top to bottom.
    """
    matrix = init_matrix(sequence1, sequence2, gap_opening_cost)
    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            matrix[i+1][j+1] = new_value_computation(sequence1[i], sequence2[j],
                                                     gap_opening_cost, substitution_cost,
                                                     matrix[i][j], matrix[i][j+1], matrix[i+1][j])
            
    return matrix

#function to compute the traceback (only one traceback is computed)
def traceback(matrix, sequence1, sequence2, gap_opening_cost,  substitution_cost):
    """
    Implement the traceback part of the algorithm
    With the given matrix traceback which cells were taken to complete the path from 
    the top left corner to the bottom right corner.
    """
    traceback = [(len(sequence1),len(sequence2))]
    i = len(sequence1) 
    j = len(sequence2) 

    while(i > 0 and j > 0):
        if matrix[i][j] == matrix[i-1][j-1] + substitution_cost[(sequence1[i-1], sequence2[j-1])]:
            traceback.append((i-1,j-1))
            #print('diag')
            i = i-1
            j = j-1
        elif matrix[i][j] == matrix[i-1][j] + gap_opening_cost:
            traceback.append((i-1,j))
            #print('top')
            i = i-1
        else:
            traceback.append((i,j-1))
            #print('left')
            j = j-1
            
    return traceback

#to create aligned sequences
def alingment_build(traceback, seq1, seq2):
    """
    Implement the alingment creation.
    Given the traceback figure out which editing operations were used to create the alingment.
    """
    seq1_align = ''
    seq2_align = ''
    j = k = 0
    
    for i in range(len(traceback)-1, 0, -1):
        if traceback[i][0]+1 == traceback[i-1][0] and traceback[i][1] +1 == traceback[i-1][1]:
            #print('diag')
            seq1_align += seq1[j]
            seq2_align += seq2[k]
            j += 1
            k += 1
        elif traceback[i][0]+1 == traceback[i-1][0] and traceback[i][1] == traceback[i-1][1]:
            #print('top')
            seq1_align += seq1[j]
            seq2_align += '_'
            j += 1
        else:
            #print('left')
            seq1_align += '_'
            seq2_align += seq2[k]
            k += 1
    
    return seq1_align, seq2_align

needleman_wunsch("data/s1.fasta", "data/s2.fasta", "data/blosum62.txt", -6)
needleman_wunsch("data/s2.fasta", "data/s3.fasta", "data/pam250.txt", -8)