
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
    
    
def visualize_matrix(matrix):
    """
    Implement the visualization of a matrix.
    Can be used for self check
    """
    for i in range(len(matrix)):
        print(*matrix[i], sep='\t')
        

def complete_d_p_q_computation(seq_1, seq_2, cost_gap_open, cost_gap_extend, substitutions=None):
    """
    Implement the recursive computation of matrices D, P and Q
    """
    matrix_d = init_matrix_d(seq_1, seq_2, cost_gap_open, cost_gap_extend)
    matrix_p = init_matrix_p(seq_1, seq_2)
    matrix_q = init_matrix_q(seq_1, seq_2)
    
    if substitutions == None:
        for i in range(1,len(seq_1)+1):
            for j in range(1,len(seq_2)+1):
                matrix_p[i][j] = max(matrix_d[i-1][j] + cost_gap_open + cost_gap_extend,
                                    matrix_p[i-1][j] + cost_gap_extend)
                matrix_q[i][j] = max(matrix_d[i][j-1] + cost_gap_open + cost_gap_extend,
                                    matrix_q[i][j-1] + cost_gap_extend)
                if seq_1[i-1] == seq_2[j-1]:
                    matrix_d[i][j] = max(matrix_d[i-1][j-1] + 1, matrix_p[i][j], matrix_q[i][j])
                else:
                    matrix_d[i][j] = max(matrix_d[i-1][j-1] - 1, matrix_p[i][j], matrix_q[i][j])
                    
    else:
        for i in range(1,len(seq_1)+1):
            for j in range(1,len(seq_2)+1):
                matrix_p[i][j] = max(matrix_d[i-1][j] + cost_gap_open + cost_gap_extend,
                                    matrix_p[i-1][j] + cost_gap_extend)
                matrix_q[i][j] = max(matrix_d[i][j-1] + cost_gap_open + cost_gap_extend,
                                    matrix_q[i][j-1] + cost_gap_extend)
                matrix_d[i][j] = max(matrix_d[i-1][j-1] + substitutions[(seq_1[i-1], seq_2[j-1])],
                                     matrix_p[i][j], matrix_q[i][j])
                
    return matrix_d, matrix_p, matrix_q


"""
You are working with 3 matrices simultaneously.
You can store your path as a list of cells.
A cell can be a tuple: coordinates, matrix_name.
And coordinates is a tuple of indexex i, j.

Cell example: ((0, 2), "d")
Path example: [((2, 4), 'd'), ((2, 4), 'q'), ((2, 3), 'q'), ((2, 2), 'd'), ((1, 1), 'd'), ((0, 0), 'd')]

"""

def compute_all_tracebacks(seq1, seq2, d_matrix, p_matrix, q_matrix,
                           cost_gap_open, cost_gap_extend, substitution=None):
    """
    Implement a search for all possible paths from the bottom right corner to the top left.
    Implement 'find_all_previous' and check_complete first.
   
    """
    
    all_paths = [[((len(seq1),len(seq2)), 'd')]]
    check_path = False
    i = 0
    
    while i < len(all_paths):
        while(check_path ==  check_complete(all_paths[i])):
            cells = find_all_previous(all_paths[i][-1], seq1, seq2, d_matrix, p_matrix, q_matrix,
                                      cost_gap_open, cost_gap_extend, substitution)
            
            if len(cells) == 1:
                all_paths[i].append(cells[0])
            else:
                for j in range(len(cells)-1):
                    all_paths.append(all_paths[i][:])
                    all_paths[-1].append(cells[j+1])
                all_paths[i].append(cells[0])
        
        i += 1
            
    return all_paths
    
def find_all_previous(cell, seq1, seq2, d_matrix, p_matrix, q_matrix,
                   cost_gap_open, cost_gap_extend):
    parent_cells = []
    """
    Implement a search for all possible previous cells.
    """
    if cell[1] == 'p':
        if p_matrix[cell[0][0]][cell[0][1]] == d_matrix[cell[0][0]-1][cell[0][1]] + cost_gap_open + cost_gap_extend:
            parent_cells.append(((cell[0][0]-1, cell[0][1]), 'd'))
        if p_matrix[cell[0][0]][cell[0][1]] == p_matrix[cell[0][0]-1][cell[0][1]] + cost_gap_extend:
            parent_cells.append(((cell[0][0]-1, cell[0][1]), 'p'))
    elif cell[1] == 'q':
        if q_matrix[cell[0][0]][cell[0][1]] == d_matrix[cell[0][0]][cell[0][1]-1] + cost_gap_open + cost_gap_extend:
            parent_cells.append(((cell[0][0], cell[0][1]-1), 'd'))
        if q_matrix[cell[0][0]][cell[0][1]] == q_matrix[cell[0][0]][cell[0][1]-1] + cost_gap_extend:
            parent_cells.append(((cell[0][0], cell[0][1]-1), 'q'))
    else:
        if cell[0][0] == 0:
            parent_cells.append(((cell[0][0], cell[0][1]-1), 'd'))
            return parent_cells
        if cell[0][1] == 0:
            parent_cells.append(((cell[0][0]-1, cell[0][1]), 'd'))
            return parent_cells
        if d_matrix[cell[0][0]][cell[0][1]] == p_matrix[cell[0][0]][cell[0][1]]:
            parent_cells.append(((cell[0][0], cell[0][1]), 'p'))
        if d_matrix[cell[0][0]][cell[0][1]] == q_matrix[cell[0][0]][cell[0][1]]:
            parent_cells.append(((cell[0][0], cell[0][1]), 'q'))
        if seq1[cell[0][0]-1] == seq2[cell[0][1]-1]:
            if d_matrix[cell[0][0]][cell[0][1]] == d_matrix[cell[0][0]-1][cell[0][1]-1] + 1:
                parent_cells.append(((cell[0][0]-1, cell[0][1]-1), 'd'))
        else:
            if d_matrix[cell[0][0]][cell[0][1]] == d_matrix[cell[0][0]-1][cell[0][1]-1] - 1:
                parent_cells.append(((cell[0][0]-1, cell[0][1]-1), 'd'))
            
    return parent_cells

def check_complete(path):
    """
    Implement a function which checks if the traceback path is complete.
    """
    if path[-1][0]== (0,0):
        return True
    else:
        return False
        
              
 
def alignment(traceback_path, seq1, seq2):
    """
    Implement creation of the alignment with given traceback path and sequences1 and 2
    """
    alignment_seq1 = ''
    alignment_seq2 = ''
    j = 0
    k = 0
    i = len(traceback_path)-2
    
    while i > -1:
        if traceback_path[i][1] == 'd' and traceback_path[i][0][0] == 0:
            alignment_seq1 += '-'
            alignment_seq2 += seq2[k]
            k += 1
            i -= 1
        elif traceback_path[i][1] == 'd' and traceback_path[i][0][1] == 0:
            alignment_seq1 += seq1[j]
            alignment_seq2 += '-'
            j += 1
            i -= 1
        elif traceback_path[i][1] == 'd':
            #print('d')
            alignment_seq1 += seq1[j]
            alignment_seq2 += seq2[k]
            j += 1
            k += 1
            i -= 1
        elif traceback_path[i][1] == 'p':
            #print('p')
            alignment_seq1 += seq1[j]
            alignment_seq2 += '-'
            j += 1
            i -= 1
            if traceback_path[i][1] == 'd':
                i -= 1
        else:
            #print('q')
            alignment_seq1 += '-'
            alignment_seq2 += seq2[k]
            k += 1
            i -= 1
            if traceback_path[i][1] == 'd':
                i -= 1
    return alignment_seq1, alignment_seq2
    
          
