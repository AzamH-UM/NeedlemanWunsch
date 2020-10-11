/** Azam Hussain
 * NW Algorithm
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>

using namespace std;

struct cell {
   int NWScore = 0;
   int Arrow;
   int gapLenTop = 1;
   int gapLenLeft = 1;
   int preTop;
   int preLeft; 
};

string readFromPdb(string pdb);
string readFromFasta(string fasta);

int main(int argc, char* argv[]) {

    string seq1;
    string seq2;
    
    //CHANGES HERE
    bool pdb = false; //change to read from pdb file

    if (argc > 1) {
        string fasta1 = argv[1];
        string fasta2 = argv[2];
        seq1 = readFromFasta(fasta1);
        seq2 = readFromFasta(fasta2);
    }

    else if (pdb) {
        seq2 = readFromPdb("./1ko3_pdb.pdb");
        seq1 = readFromPdb("./3v2b_pdb.pdb");
    }






    else {
        seq1 = "FSSSELYNWFTLTNLKPDANTGVVNFDIPGYIHDFASKDATVTLASNPLSWLVAATGWHYGEVDLCISWSRSKQAQAQEGSVSITTNYRDWGAYWQGQARIYDLRRTEAEIPIFLGSYAGATPSGALGKQNYVRISIVNAKDIVALRVCLRPKSIKFWGRSATLF";
        seq2 = "GKRILLLEKERNLAHFLSLELQKEQYRVDLVEEGQKALSMALQTDYDLILLNVNLGDMMAQDFAEKLSRTKPASVIMILDHWEDLQEELEVVQRFAVSYIYKPVLIENLVARISAIFRGRDFI";
    }
    int indel = -11; //gap opening (insertion/deletion)
    int extension = -1; //gap extension
    //END OF CHANGES

    unordered_map<string, int> AminoToInt({ //maps single letter amino acid to int for BLOSUM matrix
		{"A", 0},
		{"R", 1},
		{"N", 2},
		{"D", 3},
		{"C", 4},
		{"Q", 5},
		{"E", 6},
		{"G", 7},
		{"H", 8},
		{"I", 9},
		{"L", 10},
		{"K", 11},
		{"M", 12},
		{"F", 13},
		{"P", 14},
		{"S", 15},
		{"T", 16},
		{"W", 17},
		{"Y", 18},
		{"V", 19},
		});


    int BLOSUM62[20][20] = { //scores for pairwise alignment, amino acids in same order as map
        { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
        {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
        {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
        {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
        { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
        {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
        {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
        { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
        {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
        {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
        {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
        {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
        {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
        {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
        { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
        { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
        {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
        { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}
    };
    
    cell NWArray[seq1.length() + 1][seq2.length() + 1];
    int finalscore;

    //initiates the first cells
    NWArray[0][0].NWScore = 0;
    NWArray[0][0].preLeft = -1000;
    NWArray[0][0].preTop = -1000;


    NWArray[1][0].NWScore = indel;
    NWArray[1][0].preLeft = -1000;
    NWArray[1][0].preTop = -1000;

    NWArray[0][1].NWScore = indel;
    NWArray[0][1].preTop = -1000;
    NWArray[0][1].preLeft = -1000;
    
    
    //initiates first col to gap extension scores;
    for (int i = 2; i <= seq1.length(); i++) {
        cell *ci0 = &NWArray[i][0];
        (*ci0).NWScore = NWArray[i-1][0].NWScore+ extension;
        (*ci0).preTop = -1000;
        (*ci0).preLeft = -1000;
        (*ci0).gapLenLeft = 0;
    }

    //initates the first row to gap extension scores;
    for (int j = 2; j <= seq2.length(); j++) {
        cell *c0j = &NWArray[0][j];
        (*c0j).NWScore = NWArray[0][j-1].NWScore+ extension;
        (*c0j).preLeft = -1000;
        (*c0j).preTop = -1000;
        (*c0j).gapLenTop = 0;
    }

     
    int top; //top cell, represents indel
    int left; //bot cell, represents indel
    int match; //diagonal cell, represent match/mismatch

    //nested for loop to start checking 3 neighbors of each matrix element
    for (int i = 1; i <= seq1.length(); i++) { //i for rows
        for (int j = 1; j <= seq2.length(); j++) { //j for columns
        
            cell *cij = &NWArray[i][j]; //current cell at i,j


            // Gap checking left
            int left1 = NWArray[i][j-1].NWScore + indel;   //opening a new gap
            int left2 = NWArray[i][j-1].preLeft + extension; //extending the previous gap

            if (left1 > left2) {
                left = left1; //opening a new gap is cheaper
            }
            else {
                left = left2; //extending the previous gap is better
                (*cij).gapLenLeft = NWArray[i][j-1].gapLenLeft + 1; //add one to the length of the gap

            }

            // Gap checking top
            int top1 = NWArray[i-1][j].NWScore + indel;
            int top2 = NWArray[i-1][j].preTop + extension;

            if (top1 > top2) {
                top = top1; //opening a gap is better
            }
            else {
                top = top2; //extending the previous gap is better
                (*cij).gapLenTop = NWArray[i-1][j].gapLenTop + 1;
            }

            (*cij).preLeft = left; //store left score of this cell
            (*cij).preTop = top;   //store top score of this cell

            int res1 = AminoToInt[seq1.substr(i-1,1)]; //Find ith amino acid of seq1
            int res2 = AminoToInt[seq2.substr(j-1,1)]; //Find jth amino acid of seq2
            int matchscore = BLOSUM62[res1][res2];     //Get pairwise (ith,jth) score
            match = NWArray[i-1][j-1].NWScore + matchscore;

            if (match > top && match > left) {
                (*cij).NWScore = match;
                (*cij).Arrow = 1;
                finalscore = match;
            }
            else if (left > top) {
                (*cij).NWScore = left;
                (*cij).Arrow = 2;
                finalscore = left;
            }
            else {
                (*cij).NWScore = top;
                (*cij).Arrow = 3;
                finalscore = top;
            }
            

            
        }
    }
    
    //Alignment: Start from corner and follow path made by "arrows" to preTopious cell
    bool done = false;
    string align1 = ""; //output for sequence1
    string align2 = ""; //output for sequence2
    string col = "";    //output for colon denoting match
    int i = seq1.length(); 
    int j = seq2.length();
    int matching = 0;
    

    while(!done){

        cell *cij = &NWArray[i][j];

        if (i == 0) {           //stuck on the top most row
            align1.insert(0, "-");               
            align2.insert(0, seq2.substr(j-1,1));
            j--;
            col.insert(0, " ");
        }
        else if (j == 0) {     //stuck on the left most column
            align1.insert(0, seq1.substr(i-1,1));
            align2.insert(0, "-");
            i--;
            col.insert(0, " ");
        }

        else if ((*cij).Arrow == 1) {        //match
            align1.insert(0, seq1.substr(i-1,1));
            align2.insert(0, seq2.substr(j-1,1));

            if (seq1.substr(i-1,1).compare(seq2.substr(j-1,1)) == 0) {
                col.insert(0, ":");
                matching++;
            }
            else {
                col.insert(0, " ");
            }
            j--;
            i--;
            
        }
        
        else if ((*cij).Arrow == 2){ //left gap
            //follow gap
            int gap = (*cij).gapLenLeft;
            for (int n = 0; n < gap; n++) {
                if (j > 0) {
                    align1.insert(0, "-");               
                    align2.insert(0, seq2.substr(j-1,1));
                    j--;
                    col.insert(0, " ");
                }
            }
            
        }
        else {   //top gap               
            //follow gap
            int gap = (*cij).gapLenTop;
            for (int n = 0; n < gap; n++) {
                if (i > 0) {
                    align1.insert(0, seq1.substr(i-1,1));
                    align2.insert(0, "-");
                    i--;
                    col.insert(0, " ");
                }
            
            }
        }    
        
        if (i == 0 && j == 0) {                 //beginning and termination of loop
            done = true;
        }

        
    }

    int score = 0;

    //score calculation:
    for (int i = 0; i < align1.length(); i++) {
        string item1 = align1.substr(i,1);
        string item2 = align2.substr(i,1);

        if (item1.compare("-") == 0) {
            score += indel;
            bool gap = true;
            while (gap) {
                i++;
                item1 = align1.substr(i,1);
                if (item1.compare("-") == 0) {
                    score += extension;
                }
                else {
                    gap = false;
                    i--;
                }
            } 
        }
        else if (item2.compare ("-") == 0) {
            score += indel;
            bool gap = true;
            while (gap) {
                i++;
                item2 = align2.substr(i,1);
                if (item2.compare("-") == 0) {
                    score += extension;
                }
                else {
                    gap = false;
                    i--;
                }
            } 
        }
        else {
            score += BLOSUM62[AminoToInt[item1]][AminoToInt[item2]];
        }

    }

    /*
    for (int i = 0; i <= seq1.length(); i++) { //i for rows
        for (int j = 0; j <= seq2.length(); j++) { //j for columns
            cout << NWArray[i][j].NWScore << " ";
        }
        cout << endl;
    }
    */

    
    //cout << align1 << "\n" << col << "\n" << align2 << "\n" << "Identical Length: " << matching << "\n";
    cout << align1 << "\n" << align2 << "\n";
    cout << "AlignmentScore: " << score << "\nMatrixScore: " << finalscore << "\n";
    
    return 0;
}


string readFromPdb(string pdb) {
    //AMINO ACIDS
	unordered_map<string, string> AminoMap({ 
		{"ALA","A"},
		{"ARG","R"},
		{"ASN","N"},
		{"ASP","D"},
		{"CYS","C"},
		{"GLU","E"},
		{"GLN","Q"},
		{"GLY","G"},
		{"HIS","H"},
		{"ILE","I"},
		{"LEU","L"},
		{"LYS","K"},
		{"MET","M"},
		{"PHE","F"},
		{"PRO","P"},
		{"SER","S"},
		{"THR","T"},
		{"TRP","W"},
		{"TYR","Y"},
		{"VAL","V"},
		});

	string line; //holding new lines of pdb file for getline()
	ifstream pdbReader(pdb);
	
	string currentRes;
	string checkedRes;
    string sequence = "";
	
	if (pdbReader.is_open()) {
		while (getline(pdbReader, line)) {
			currentRes = line.substr(23, 3);
			if (line.substr(0, 4).compare("ATOM") == 0 && currentRes.compare(checkedRes) != 0) { //Checks PDB for ATOM containing amino acid sequence
				sequence.append(AminoMap[line.substr(17,3)]);          //Written to fasta file
				checkedRes = currentRes;
			}
		}
		pdbReader.close();
	}
	else cout << "Unable to open file";
    return sequence;

}
 
string readFromFasta(string fasta) {
    string line; //holding new lines of pdb file for getline()
	ifstream fastaReader(fasta);
    getline(fastaReader, line);
    getline(fastaReader,line);
    return line;
}