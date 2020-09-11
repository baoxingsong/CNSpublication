package me.songbx.mains;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.KmerReadImpl;

public class CnsSeqOfGenomeKmerSpectrum {
	public static void main ( String argv[] ) {
		new CnsSeqOfGenomeKmerSpectrum( "./nonIntronUtr_OpOrTFBS.fa", "./b73_m20_mer_counts_dumps.fa", "./nonIntronUtr_OpOrTFBS.kmer_spectrum" );
		new CnsSeqOfGenomeKmerSpectrum( "./nonOpTfbsIntronUtrCns.fa", "./b73_m20_mer_counts_dumps.fa", "./nonOpTfbsIntronUtrCns.kmer_spectrum" );
	}
	public CnsSeqOfGenomeKmerSpectrum( String fastaFile, String kmerFile, String outptuFile ) {
		int k = 20;
		HashMap<String, String> sequences = new ChromoSomeReadImpl(fastaFile).getChromoSomeHashMap();
		HashMap<String, Integer> kMers = new HashMap<String, Integer>();
		for( String name : sequences.keySet() ) {
			String seq = sequences.get(name);
			for ( int i=0; i<= seq.length()-k; ++i ) {
				String mer = seq.substring(i, i+k);
				String rMer = ChromoSomeReadImpl.getReversecomplementary(mer);
				if( rMer.compareTo(mer) < 0 ) {
					mer = rMer;
				}
				if( kMers.containsKey(mer) ) {
					kMers.put(mer, kMers.get(mer)+1);
				}else {
					kMers.put(mer, 1);
				}
			}
		}
		sequences.clear();
		HashMap<String, Integer> genomeKmers = (new KmerReadImpl(kmerFile)).getKmers();
		try {
			PrintWriter outPut2 = new PrintWriter(outptuFile);
			outPut2.println("kmer\tkmerFrequency\tgenomeKmerFrequency\tkmerGcContent");
			for( String kmer : kMers.keySet() ) {
				outPut2.println(kmer + "\t" + kMers.get(kmer) + '\t' + genomeKmers.get(kmer) + "\t" +KmerReadImpl.getGcContent(kmer));
			}
			outPut2.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
}
