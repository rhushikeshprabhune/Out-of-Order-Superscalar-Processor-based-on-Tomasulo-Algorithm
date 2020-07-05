#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
using namespace std;

enum {
  IF = 0,
  ID,
  IS,
  EX,
  WB
};

typedef struct  
{
    unsigned int t_tag, t_LRU_counter, t_dirty_bit, t_valid_bit; //Define structure for 
}table;

struct fake_rob{
    signed int optype, stage, tag, src1, src2, dest, src1_rdy, src2_rdy, src1_name, src2_name, addr, cycle_counter, EX_latency;
    int startIF,endIF,startID,endID,startIS,endIS,startEX,endEX,startWB, endWB;
};

struct rmt{
    int reg_name, reg_rdy;
};

class L2CACHE
{
	public:
		unsigned int BLOCKSIZE, L2_SIZE, L2_ASSOC, sets, bo, indexsizeL2, tagsizeL2, setsL2;
	int L2_latency;
	unsigned int readsL2, writesL2, read_missesL2, write_missesL2, write_backsL2, sector_miss, block_miss;
	table ** ArrL2; //Two dimentional array to hold the index values which will point to the blocks in that index


	L2CACHE(unsigned int c_BLOCKSIZE, unsigned int c_L2_SIZE, unsigned int c_L2_ASSOC)
	{

		BLOCKSIZE = c_BLOCKSIZE;
		L2_SIZE = c_L2_SIZE;
		L2_ASSOC = c_L2_ASSOC;
		if (c_L2_SIZE > 0)
		{
			sets = c_L2_SIZE / (c_BLOCKSIZE * c_L2_ASSOC); //number of sets
			bo = log2(c_BLOCKSIZE); //block offset
			indexsizeL2 = log2(sets); //number of bits in index value
			tagsizeL2 = 32 - bo - indexsizeL2;
			readsL2 = 0;
			writesL2 = 0;
			read_missesL2 = 0;
			write_missesL2 = 0;
			write_backsL2 = 0;
			sector_miss = 0;
			block_miss = 0;
			L2_latency = 0;
			ArrL2 = new table * [sets]; //Define the two dimentional array
			for (unsigned int i = 0; i < sets; i++)
				ArrL2[i] = new table[L2_ASSOC];

			for (unsigned int i = 0; i < sets; i++)
			{
				for (unsigned int j = 0; j < L2_ASSOC; j++)
				{
					ArrL2[i][j].t_LRU_counter = j; //initialise the lru counters, dirty bits,valid bits and tags
					ArrL2[i][j].t_dirty_bit = 0;
					ArrL2[i][j].t_valid_bit = 0;
					ArrL2[i][j].t_tag = 0;
				}
			}
		}
	}


	int readFromAddressL2(unsigned int);
	void writeToAddressL2(unsigned int);
	void printStatsL2();
};

class L1CACHE
{
	public:
		unsigned int BLOCKSIZE, L1_SIZE, L1_ASSOC, sets, block, bo, indexsize, tagsize;
	int L1_latency;
	unsigned int reads, writes, read_misses, write_misses, write_backs;
	L2CACHE * next; //Next level cache
	table ** Arr; //Two dimentional array to hold the index values which will point to the blocks in that index
	int l2count;

	L1CACHE(unsigned int c_BLOCKSIZE, unsigned int c_L1_SIZE, unsigned int c_L1_ASSOC, L2CACHE * c_next)
	{
		BLOCKSIZE = c_BLOCKSIZE;
		L1_SIZE = c_L1_SIZE;
		L1_ASSOC = c_L1_ASSOC;
		if (c_L1_SIZE != 0)
		{
			block = c_L1_SIZE / c_BLOCKSIZE; //number of blocks
			sets = c_L1_SIZE / (c_BLOCKSIZE * c_L1_ASSOC); //number of sets
			bo = log2(c_BLOCKSIZE); //block offset
			indexsize = log2(sets); //number of bits in index value
			tagsize = 32 - bo - indexsize; //tag length
			next = c_next; //Next level(NULL for Project:1A)                                
			reads = 0;
			writes = 0;
			read_misses = 0;
			write_misses = 0;
			write_backs = 0;
			L1_latency = 0;
			Arr = new table * [sets]; //Define the two dimentional array
			for (unsigned int i = 0; i < sets; i++)
				Arr[i] = new table[L1_ASSOC];

			for (unsigned int i = 0; i < sets; i++)
			{
				for (unsigned int j = 0; j < L1_ASSOC; j++)
				{
					Arr[i][j].t_LRU_counter = j; //initialise the lru counters, dirty bits,valid bits and tags
					Arr[i][j].t_dirty_bit = 0;
					Arr[i][j].t_valid_bit = 0;
					Arr[i][j].t_tag = 0;
				}
			}
		}
	}
	int readFromAddressL1(unsigned int);
	void writeToAddressL1(unsigned int);
	void printStatsL1();
};

int L2CACHE::readFromAddressL2(unsigned int addr)
{
	readsL2++;
	int r_hit_block = 0;
	unsigned int tag = addr >> (bo + indexsizeL2); //extract the tag from input address
	unsigned int index = ((1 << indexsizeL2) - 1) & (addr >> (bo));

	int r_hit = 0;
	for (unsigned int j = 0; j < L2_ASSOC; j++)
	{
		if (ArrL2[index][j].t_tag == tag) //find if hit or miss and if hit find the hit block
		{
			r_hit = 1;
			r_hit_block = j;
		}
	}

	if (r_hit == 1) //hit case
	{
		L2_latency = 1;
		unsigned int oldLRU;
		oldLRU = ArrL2[index][r_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L2_ASSOC; i++)
		{
			if (ArrL2[index][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
				ArrL2[index][i].t_LRU_counter += 1; //is less than hit block counter value
		}
		ArrL2[index][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
	}


	if (r_hit == 0) //miss case
	{
		L2_latency = 0;
		read_missesL2 = read_missesL2 + 1;

		int invalid = 0;
		for (unsigned int vl = 0; vl < L2_ASSOC; vl++)
		{
			if (ArrL2[index][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;
			}
		}

		if (invalid == 1) //case when only one block is empty
		{
			for (unsigned int i = 0; i < L2_ASSOC; i++)
			{
				if (ArrL2[index][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit 
				{
					for (unsigned int j = 0; j < L2_ASSOC; j++)
					{
						ArrL2[index][j].t_LRU_counter++;
					}

					ArrL2[index][i].t_tag = tag;
					ArrL2[index][i].t_LRU_counter = 0;
					ArrL2[index][i].t_valid_bit = 1;
				}
			}

		}

		if (invalid == 0) //when all are filled, check highest LRU_count and update tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = ArrL2[index][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L2_ASSOC; j++)
			{
				if (ArrL2[index][j].t_LRU_counter > temp)
				{
					temp = ArrL2[index][j].t_LRU_counter;
					highest_LRU = j;
				}
			}

			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{

				ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag;
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			if (ArrL2[index][highest_LRU].t_dirty_bit == 1)
				write_backsL2++;
			ArrL2[index][highest_LRU].t_dirty_bit = 0;

		}

		if (invalid != 0 && invalid != 1) //case when multiple blocks are empty
		{
			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_valid_bit == 0 && ArrL2[index][k].t_LRU_counter > currentHighest) //check which are empty and find the
				{ //one with highest lru count
					currentHighest = ArrL2[index][k].t_LRU_counter;
					highest_LRU = k;
				}

			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_LRU_counter < ArrL2[index][highest_LRU].t_LRU_counter)
					ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag; //update tag,lru counter and valid bit 
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			ArrL2[index][highest_LRU].t_valid_bit = 1;
		}

	}
	return (L2_latency);
}

int L1CACHE::readFromAddressL1(unsigned int addr)
{
	reads++;
	int r_hit_block = 0;
	unsigned int tagr = addr >> (bo + indexsize); //extract the tag from input address
	unsigned int indexr = ((1 << indexsize) - 1) & (addr >> (bo));
	unsigned int l2HitOrMiss = 0;

	int r_hit = 0;
	for (unsigned int j = 0; j < L1_ASSOC; j++)
	{
		if (Arr[indexr][j].t_tag == tagr) //find if hit or miss and if hit find the hit block
		{
			r_hit = 1;
			r_hit_block = j;
		}
	}

	if (r_hit == 1) //hit case
	{
		L1_latency = 5;
		unsigned int oldLRU;
		oldLRU = Arr[indexr][r_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L1_ASSOC; i++)
		{
			if (Arr[indexr][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
				Arr[indexr][i].t_LRU_counter += 1; //is less than hit block counter value
		}
		Arr[indexr][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
	}


	if (r_hit == 0) //miss case
	{
		//l2HitOrMiss = next->readFromAddressL2(addr);

		read_misses++;
		int invalid = 0;
		for (unsigned int vl = 0; vl < L1_ASSOC; vl++)
		{
			if (Arr[indexr][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;

			}
		}


		if (invalid == 1) //case when only one block is empty
		{
			for (unsigned int i = 0; i < L1_ASSOC; i++)
			{
				if (Arr[indexr][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit 
				{
					for (unsigned int j = 0; j < L1_ASSOC; j++)
					{

						Arr[indexr][j].t_LRU_counter++;
					}

					Arr[indexr][i].t_tag = tagr;
					Arr[indexr][i].t_LRU_counter = 0;
					Arr[indexr][i].t_valid_bit = 1;
				}
			}
			if (next - > L2_SIZE != 0)
				l2HitOrMiss = next - > readFromAddressL2(addr);
		}

		if (invalid == 0) //when all are filled, check highest LRU_count and update tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = Arr[indexr][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L1_ASSOC; j++)
			{
				if (Arr[indexr][j].t_LRU_counter > temp)
				{
					temp = Arr[indexr][j].t_LRU_counter;
					highest_LRU = j;
				}
			}
			unsigned int z = (Arr[indexr][highest_LRU].t_tag << indexsize | indexr) << bo;


			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{

				Arr[indexr][k].t_LRU_counter++;
			}

			Arr[indexr][highest_LRU].t_tag = tagr;
			Arr[indexr][highest_LRU].t_LRU_counter = 0;
			if (Arr[indexr][highest_LRU].t_dirty_bit == 1)
			{
				write_backs++;
				//next->writeToAddressL2(z);

			}
			Arr[indexr][highest_LRU].t_dirty_bit = 0;

			if (next - > L2_SIZE != 0)
				l2HitOrMiss = next - > readFromAddressL2(addr);
		}

		if (invalid != 0 && invalid != 1) //case when multiple blocks are empty
		{

			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[indexr][k].t_valid_bit == 0 && Arr[indexr][k].t_LRU_counter > currentHighest) //check which are empty and find the
				{ //one with highest lru count
					currentHighest = Arr[indexr][k].t_LRU_counter;
					highest_LRU = k;
				}

			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[indexr][k].t_LRU_counter < Arr[indexr][highest_LRU].t_LRU_counter)
					Arr[indexr][k].t_LRU_counter++;
			}

			Arr[indexr][highest_LRU].t_tag = tagr; //update tag,lru counter and valid bit 
			Arr[indexr][highest_LRU].t_LRU_counter = 0;
			Arr[indexr][highest_LRU].t_valid_bit = 1;

			if (next - > L2_SIZE != 0)
				l2HitOrMiss = next - > readFromAddressL2(addr);

		}

		if (next - > L2_SIZE == 0)
			L1_latency = 20;
		else if (next - > L2_SIZE > 0)
		{
			if (l2HitOrMiss == 1)
				L1_latency = 10;
			else if (l2HitOrMiss == 0)
				L1_latency = 20;
		}
	}

	return (L1_latency);
}
void L1CACHE::printStatsL1()
{
	cout << "L1 CACHE CONTENTS" << endl;
	cout << "a. number of accesses :" << reads << endl;
	cout << "b. number of misses :" << read_misses << endl;
	for (unsigned int i = 0; i < sets; i++)
	{
		cout << "set  " << dec << i << ": ";
		for (unsigned int j = 0; j < L1_ASSOC; j++)
		{
			cout << " " << hex << Arr[i][L1_ASSOC - j - 1].t_tag;
		}
		cout << '\n';
	}
}

void L2CACHE::printStatsL2()
{
	cout << "L2 CACHE CONTENTS" << endl;
	cout << "a. number of accesses :" << dec << readsL2 << endl;
	cout << "b. number of misses :" << dec << read_missesL2 << endl;
	for (unsigned int i = 0; i < sets; i++)
	{
		cout << "set  " << dec << i << ": ";
		for (unsigned int j = 0; j < L2_ASSOC; j++)
		{
			cout << " " << hex << ArrL2[i][L2_ASSOC - j - 1].t_tag;
		}
		cout << '\n';
	}
	cout << "" << endl;
}

unsigned int addr;
int cycle = 0;
int trace_fully_read = 0;
int dispatch_counter = 0;
int schedule_counter = 0;
int instr_count = 0;
int rob_counter = 0;
fake_rob ROB[1024];
int head = 0;
int tail = 0;
rmt RMT[128];
ifstream infile;

void Fetch(int N)
{
	signed int optype, dest, src1, src2, addr, PC;
	int diff = 0;
	diff = 2 * N - dispatch_counter;
	int count = 0;
	//if dispatch queue is filled more than N, fill till it becomes 2N
	if (diff < N)
	{
		for (int i = 0; i < diff; i++)
		{
			if (infile >> hex >> PC >> dec >> optype >> dest >> src1 >> src2 >> hex >> addr >> dec)
			{
				ROB[tail].addr = addr;
				// cout << addr << endl; 
				ROB[tail].optype = optype;
				ROB[tail].src1 = src1;
				ROB[tail].src2 = src2;
				ROB[tail].dest = dest;
				ROB[tail].stage = IF;
				ROB[tail].startIF = cycle;
				ROB[tail].tag = instr_count;
				tail++;
				if (tail == 1024)
					tail = 0;
				dispatch_counter++;
				instr_count++;
				rob_counter++;
			}
			else
			{
				trace_fully_read = 1;
			}
		}
	}
	//if dispatch queue is filled less than N, can put max BW number of instrs i.e N
	else
	{
		for (int i = 0; i < N; i++)
		{
			if (infile >> hex >> PC >> dec >> optype >> dest >> src1 >> src2 >> hex >> addr >> dec)
			{
				ROB[tail].addr = addr;
				ROB[tail].optype = optype;
				ROB[tail].src1 = src1;
				ROB[tail].src2 = src2;
				ROB[tail].dest = dest;
				ROB[tail].stage = IF;
				ROB[tail].startIF = cycle;
				ROB[tail].tag = instr_count;
				tail++;
				if (tail == 1024)
					tail = 0;
				dispatch_counter++;
				instr_count++;
				rob_counter++;

			}
			else
			{
				trace_fully_read = 1;
			}
		}
	}
}

bool AdvanceCycle()
{
	if (trace_fully_read == 1 && rob_counter == 0)
	{
		return false;
	}
	else
	{
		cycle++;
		return true;
	}
}

void Dispatch(int S, int N)
{

	int i = head;
	int num_dispatch = 0;
	int diff_sch = S - schedule_counter;

	if (diff_sch != 0)
	{
		while (i != tail && num_dispatch < diff_sch)
		{
			if (ROB[i].stage == ID)
			{

				ROB[i].stage = IS;
				ROB[i].endID = cycle;
				ROB[i].startIS = cycle;

				if (ROB[i].src1 != -1)
				{
					ROB[i].src1_rdy = RMT[ROB[i].src1].reg_rdy;
				}
				else if (ROB[i].src1 == -1)
				{
					ROB[i].src1_rdy = 1;
				}

				if (ROB[i].src2 != -1)
				{
					ROB[i].src2_rdy = RMT[ROB[i].src2].reg_rdy;
				}
				else if (ROB[i].src2 == -1)
				{
					ROB[i].src2_rdy = 1;
				}

				if (ROB[i].src1 != -1)
					ROB[i].src1_name = RMT[ROB[i].src1].reg_name;

				if (ROB[i].src2 != -1)
					ROB[i].src2_name = RMT[ROB[i].src2].reg_name;

				if (ROB[i].dest != -1)
				{
					RMT[ROB[i].dest].reg_name = ROB[i].tag;
					RMT[ROB[i].dest].reg_rdy = 0;
				}
				schedule_counter++;
				dispatch_counter--;
				num_dispatch++;
			}
			if (i != 1023)
				i++;
			else if (i == 1023)
				i = 0;
		}
	}

	int j = head;
	while (j != tail)
	{
		if (ROB[j].stage == IF)
		{
			ROB[j].stage = ID;
			ROB[j].endIF = cycle;
			ROB[j].startID = cycle;
		}
		if (j != 1023)
			j++;
		else if (j == 1023)
			j = 0;
	}
}

void Issue(L1CACHE * cache, int N)
{
	int execution_counter = 0;
	int i = head;
	while (i != tail && execution_counter < N)
	{
		if (ROB[i].stage == IS)
		{

			if (ROB[i].src1_rdy == 1 && ROB[i].src2_rdy == 1)
			{
				execution_counter++;
				schedule_counter--;
				ROB[i].stage = EX;
				ROB[i].endIS = cycle;
				ROB[i].startEX = cycle;
				if (ROB[i].optype == 0)
					ROB[i].EX_latency = 1;
				else if (ROB[i].optype == 1)
					ROB[i].EX_latency = 2;
				else if (ROB[i].optype == 2) //need to alter this optype to incorporate cache
					if (cache - > L1_SIZE != 0)
					{
						ROB[i].EX_latency = cache - > readFromAddressL1(ROB[i].addr);
					}
				else
				{
					ROB[i].EX_latency = 5;
				}
			}
		}
		i++;
		if (i > 1023)
		{
			i = 0;
		}
	}
}

void Execute()
{
	int i = head;
	while (i != tail)
	{

		if (ROB[i].stage == EX && ROB[i].EX_latency == 1)
		{
			ROB[i].stage = WB;
			ROB[i].endEX = cycle;
			ROB[i].startWB = cycle;

			if (ROB[i].tag == RMT[ROB[i].dest].reg_name)
			{
				RMT[ROB[i].dest].reg_rdy = 1;
			}

			//wake up
			int j = head;
			while (j != tail)
			{
				if (ROB[i].tag == ROB[j].src1_name)
				{
					ROB[j].src1_rdy = 1;
				}
				if (ROB[i].tag == ROB[j].src2_name)
				{
					ROB[j].src2_rdy = 1;
				}

				j++;
				if (j > 1023)
				{
					j = 0;
				}

			}
		}

		if (ROB[i].stage == EX)
		{
			ROB[i].EX_latency--;
		}

		i++;
		if (i > 1023)
		{
			i = 0;
		}

	}

}

void FakeRetire()
{

	while (head != tail && ROB[head].stage == WB)
	{
		if (ROB[head].stage == WB)
		{
			ROB[head].endWB = ROB[head].startWB + 1;
		}
		rob_counter--;
		cout << ROB[head].tag << " " << "fu{" << ROB[head].optype << "} src{" << ROB[head].src1 << "," << ROB[head].src2 << "} dst{" << ROB[head].dest << "} IF{" << dec << ROB[head].startIF << "," << ROB[head].endIF - ROB[head].startIF << "} ID{" << ROB[head].startID << "," << ROB[head].endID - ROB[head].startID << "} IS{" << ROB[head].startIS << "," << ROB[head].endIS - ROB[head].startIS << "} EX{" << ROB[head].startEX << "," << ROB[head].endEX - ROB[head].startEX << "} WB{" << ROB[head].startWB << ",1}" << endl;

		head++;
		if (head == 1024)
		{
			head = 0;
		}
	}

}

int main(int argc, char * argv[])
{

	infile.open(argv[8]);
	for (int i = 0; i < 128; i++)
	{
		RMT[i].reg_rdy = 1;
	}
	char * trace_file;
	int S, N;
	S = atoi(argv[1]);
	N = atoi(argv[2]);
	unsigned int c_BLOCKSIZE = atoi(argv[3]);
	unsigned int c_L1_SIZE = atoi(argv[4]);
	unsigned int c_L1_ASSOC = atoi(argv[5]);
	unsigned int c_L2_SIZE = atoi(argv[6]);
	unsigned int c_L2_ASSOC = atoi(argv[7]);

	L2CACHE L2(c_BLOCKSIZE, c_L2_SIZE, c_L2_ASSOC);
	L2CACHE * c_next = & L2;
	L1CACHE L1(c_BLOCKSIZE, c_L1_SIZE, c_L1_ASSOC, c_next);

	do {
		FakeRetire();
		Execute();
		Issue( & L1, N);
		Dispatch(S, N);
		Fetch(N);
	}
	while (AdvanceCycle());
	if (c_L1_SIZE != 0)
	{
		L1.printStatsL1();
		cout << '\n';
	}
	if (c_L2_SIZE != 0)
	{
		L2.printStatsL2();
	}
	printf("CONFIGURATION\n");
	printf("superscalar bandwidth (N) = %d\n", N);
	printf("dispatch queue size (2*N) = %d\n", 2 * N);
	printf("schedule queue size (S) = %d\n", S);

	printf("RESULTS\n");
	printf("number of instructions = %d\n", instr_count);
	printf("number of cycles = %d\n", cycle);
	printf("IPC = %.2f\n", ((float) instr_count) / ((float) cycle));

}
