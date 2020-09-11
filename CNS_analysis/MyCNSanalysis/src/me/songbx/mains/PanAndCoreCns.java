package me.songbx.mains;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;
import me.songbx.impl.ChromoSomeReadImpl;

/*
 * this script was designed to identify pan and core CNS
 * CNS is defined as cigar string M positions
 * 
 * **/

public class PanAndCoreCns {
	public static void main(String[] args) {	
		HashMap<String, String> genome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		ArrayList<String> bwFiles = new ArrayList<String> ();
		bwFiles.add("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5.bw");
		bwFiles.add("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5.bw");
		bwFiles.add("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw");
		bwFiles.add("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5.bw");
		bwFiles.add("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5.bw");
		for ( int i=1; i<=bwFiles.size(); ++i ) {
			ArrayList<ArrayList<Integer>> indexss = new ArrayList<ArrayList<Integer>>();
			for ( int j=0; j<bwFiles.size(); ++j ) {
				indexss.add(new ArrayList<Integer>());
				indexss.get(indexss.size()-1).add(j);
			}
			
			ArrayList<ArrayList<Integer>> newIndexss = copyt2DArrayList(indexss);
			while( indexss.get(0).size() < i ) {
				newIndexss = new ArrayList<ArrayList<Integer>>();
				for ( int j=0; j<bwFiles.size(); ++j ) {
					for( int k=0; k<indexss.size(); ++k ) {
						if ( indexss.get(k).get(indexss.get(k).size()-1) <j  ) {
							newIndexss.add(new ArrayList<Integer>());
							for( int l : indexss.get(k) ) {
								newIndexss.get(newIndexss.size()-1).add(l);
							}
							newIndexss.get(newIndexss.size()-1).add(j);
						}
					}
				}
				indexss = new ArrayList<ArrayList<Integer>>();
				indexss=copyt2DArrayList(newIndexss);
			}
			System.out.println("i=" + i + ":");
			for( int j=0; j<newIndexss.size(); ++j ) {
				ArrayList<String> wantedBwFiles = new ArrayList<String>();
				wantedBwFiles.clear();
				for( int k=0; k<newIndexss.get(j).size(); ++k ) {
					wantedBwFiles.add(bwFiles.get(newIndexss.get(j).get(k)));
				}
				System.out.print(wantedBwFiles.size());
				for( String bwf : wantedBwFiles ) {
					System.out.print(" " + bwf);
				}
				System.out.println();
				ArrayList<Integer> count = count( genome, wantedBwFiles );
				System.out.println("\t" + count.get(0) + "\t" + count.get(1)); // core CNS & pan CNS
			}
		}
	}
	public static ArrayList<ArrayList<Integer>> copyt2DArrayList( ArrayList<ArrayList<Integer>> indexss ){
		ArrayList<ArrayList<Integer>> newIndexss = new ArrayList<ArrayList<Integer>>();
		for(int i=0; i<indexss.size(); ++i ) {
			newIndexss.add(new ArrayList<Integer>() );
			for(int j=0; j<indexss.get(i).size(); ++j ) {
				newIndexss.get(newIndexss.size()-1).add(indexss.get(i).get(j));
			}
		}
		return newIndexss;
	}
	
	public static ArrayList<Integer> count( HashMap<String, String> genome, ArrayList<String> bwFiles ) {
		int totalAllCount = 0;
		int totalEverCount = 0;
		try {
			ArrayList<WigFileReader> wigs = new ArrayList<WigFileReader>();
			for ( String bwFilePath : bwFiles ) {
				Path bwFile = Paths.get(bwFilePath);
				WigFileReader wig = WigFileReader.autodetect(bwFile);
				wigs.add(wig);
			}
			for( String chr : genome.keySet()) {
				//System.out.println("chr:" + chr);
				for( int position=1; position<=genome.get(chr).length(); position+=10000 ) {
					
					ArrayList<Contig> results = new ArrayList<Contig>();
					for( WigFileReader wig : wigs ) {
						Contig result = wig.query(chr, position, position+10000);
						results.add(result);
					}
					
					for( int pi=position; pi<position+10000 && pi <= genome.get(chr).length(); ++pi ) {
						boolean allCoveraged = true;
						boolean everCoveraged = false;
						for( Contig result : results ) {
							Double thisMean = (double) result.get(pi);
							if ( thisMean.isNaN() || thisMean <= 0){
								allCoveraged=false;
							}else {
								//System.out.println(chr+":"+pi+"depth:"+thisMean);
								everCoveraged = true;
							}
						}
						if(allCoveraged) {
							totalAllCount = totalAllCount + 1;
						}
						if(everCoveraged) {
							totalEverCount = totalEverCount + 1;
						}
					}
				}
			}
			for( WigFileReader wig : wigs ) {
				wig.close();
			}
		}catch (IOException e) {
			e.printStackTrace();
		} catch (WigFileException e) {
			e.printStackTrace();
		}
		ArrayList<Integer> t = new ArrayList<Integer>();
		t.add(totalAllCount);
		t.add(totalEverCount);
		return t;
	}	
}
