from test_utils import generate_seq, mutate_seq
from seqalign import global_align, semiglobal_align, local_align
# from src import global_align0, semiglobal_align0, local_align0
import time

SEQ_LENGTH = 20
RUNS = 1

SEQ1 = generate_seq(SEQ_LENGTH)
SEQ2 = mutate_seq(SEQ1,.3)

LEN1 = len(SEQ1)
LEN2 = len(SEQ2)
GAP = -1
MATCH = 1
MISMATCH = -1

def timeTest():
	global SEQ1, SEQ2, GAP, MATCH, MISMATCH
	t0 = time.time()
	global_align(SEQ1, SEQ2, GAP, MATCH, MISMATCH)
	semiglobal_align(SEQ1, SEQ2, GAP, MATCH, MISMATCH)
	local_align(SEQ1, SEQ2, GAP, MATCH, MISMATCH)
	return time.time()-t0

# def timeTest0():
# 	global SEQ1, SEQ2, GAP, MATCH, MISMATCH
# 	t0 = time.time()
# 	global_align0(SEQ1, LEN1, SEQ2, LEN2, GAP, MATCH, MISMATCH)
# 	semiglobal_align0(SEQ1, LEN1, SEQ2, LEN2, GAP, MATCH, MISMATCH)
# 	local_align0(SEQ1, LEN1, SEQ2, LEN2, GAP, MATCH, MISMATCH)
# 	return time.time()-t0

def timeData(runs):
	total=timeTest()
	min_t = total
	max_t = total
	for _ in range(1,runs):
		current = timeTest()
		total+=current
		max_t = max(max_t, current)
		min_t = min(min_t,current)
	return total/runs,max_t,min_t


avg_time,max_time,min_time = timeData(RUNS)
print("Data over {} runs:\n\tAVG:\t{}\n\tMAX:\t{}\n\tMIN:\t{}\n".format(RUNS,avg_time,max_time,min_time))
