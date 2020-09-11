package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

// had problem with the following analysis, so did not use this one
// I used a python script shared by Katherine Mejia-Guerra

public class KmerClassification {
	public static void main( String args[] ) {
		new KmerClassification(2);
		new KmerClassification(3);
		new KmerClassification(4);
		new KmerClassification(5);
		new KmerClassification(6);
		new KmerClassification(7);
		new KmerClassification(9);
	} 
	public KmerClassification(int kmerSize) {
		//System.out.println(this.getReversecomplementary("ACGTCG"));
		HashMap<String, String> tokenIndexMap = createNewTokenIndexMap(kmerSize);
//		HashMap<String, Integer> sentences = seqToNewtokens( "cgagtctcacccaagagggttcaac", tokenIndexMap, 4);
//		for( String kmer : sentences.keySet() ) {
//			System.out.println(kmer + '\t' + sentences.get(kmer));
//		}
		HashSet<String> setenceHash = new HashSet<String>();
		for( String kmer : tokenIndexMap.keySet() ) {
			setenceHash.add(tokenIndexMap.get(kmer));
		}
		ArrayList<String> setenceList = new ArrayList<String>();
		for( String setence : setenceHash ) {
			setenceList.add(setence);
		}
		try {
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/kmerClusterCns/"+kmerSize+"merMatrix");
			outPut.print("class");
			for ( String setence : setenceList ) {
				outPut.print("\t" + setence);
			}
			outPut.println();
			File file = new File( "/media/bs674/1_8t/AndCns/kmerClusterCns/sequences_with_tag" );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				outPut.print(arrOfStr[0]);
				HashMap<String, Integer> sentences = seqToNewtokens( arrOfStr[1], tokenIndexMap, kmerSize);
				for ( String setence : setenceList ) {
					if( sentences.containsKey(setence) ) {
						outPut.print("\t" + sentences.get(setence));
					}else {
						outPut.print("\t0");
					}
				}
				outPut.println();
			}
			reader.close();
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * @param DNA
	 *            /RNA sequence
	 * @return reverse complementary DNA sequence(only T, no U)
	 */
	private String getReversecomplementary(String sequence) {
		StringBuffer reversecomplementary = new StringBuffer();
		for (int i = sequence.length() - 1; i >= 0; i--) {
			char c = sequence.charAt(i);
			if ('A' == c) {
				c = 'T';
			} else if ('T' == c) {
				c = 'A';
			} else if ('U' == c) {
				c = 'A';
			} else if ('C' == c) {
				c = 'G';
			} else if ('G' == c) {
				c = 'C';
			} else if ('R' == c) {
				c = 'Y';
			} else if ('Y' == c) {
				c = 'R';
			} else if ('K' == c) {
				c = 'M';
			} else if ('M' == c) {
				c = 'K';
			} else if ('B' == c) {
				c = 'V';
			} else if ('V' == c) {
				c = 'B';
			} else if ('D' == c) {
				c = 'H';
			} else if ('H' == c) {
				c = 'D';
			}
			reversecomplementary.append(c);
		}
		return reversecomplementary.toString();
	}
	public HashMap<String, String> createNewTokenIndexMap(int kmersize){
		HashMap<String, String> kmerIndexMap = new HashMap<String, String>();
		ArrayList<String> nucleotides = new ArrayList<String>();
		nucleotides.add("A");
		nucleotides.add("C");
		nucleotides.add("G");
		nucleotides.add("T");
		HashSet<String> kmers = new HashSet<String>();
		int length=0;
		while( length < kmersize ) {
			if( 0 == length ) {
				for( String kmer :  nucleotides ) {
					kmers.add(kmer);
				}
			}else {
				HashSet<String> newKmers = new HashSet<String>();
				for( String kmer : kmers  ) {
					for( String n :  nucleotides ) {
						newKmers.add(kmer+n);
					}
				}
				kmers = newKmers;
			}
			++length;
		}
		for( String kmer : kmers  ) {
			String rcKmer = getReversecomplementary(kmer);
			if( kmer.compareTo(rcKmer) <= 0 ) {
				kmerIndexMap.put(kmer, kmer+"N"+rcKmer);
			}else {
				kmerIndexMap.put(kmer, rcKmer+"N"+kmer);
			}
		}
		return kmerIndexMap;
	}
	
	public HashMap<String, Integer> seqToNewtokens( String seq, HashMap<String, String> kmerIndexMap, int kmersize) {
		seq = seq.toUpperCase();
		HashMap<String, Integer> sentences = new HashMap<String, Integer>();
    	if (seq.length() > 0 ) {
	    	for( int i=0; (seq.length()-i)>=kmersize ; ++i  ) {
	    		String kmer = seq.substring(i, i+kmersize);
	    		if( kmerIndexMap.containsKey(kmer) ) {
	    			if( sentences.containsKey(kmerIndexMap.get(kmer)) ) {
	    				sentences.put(kmerIndexMap.get(kmer), sentences.get(kmerIndexMap.get(kmer))+1);
	    			}else {
	    				sentences.put(kmerIndexMap.get(kmer), 1);
	    			}
	    		}
	    	}
	    }
	    return sentences;
	}
	
}
