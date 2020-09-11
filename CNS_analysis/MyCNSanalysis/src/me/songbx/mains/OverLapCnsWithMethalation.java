package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.model.GeneSimple;
import me.songbx.service.IfIntron;

public class OverLapCnsWithMethalation {
	public static void main(String[] args) {
		// only use those CNS never overlap genetic region
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3");
//		new OverLapCnsWithMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.wig", ifIntron);
		
		
		new OverLapCnsWithMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.wig");
		new OverLapCnsWithMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/loop.wig");
		new OverLapCnsWithMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.wig");
	}
	
	public OverLapCnsWithMethalation(String nonIntronUtrCnsWigFile, String nonOpTfbsEqtlIntronUtrCnsWigFile, IfIntron ifIntron) {
		
		int genomeLength = 0;
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		for( int i=1; i<=10; ++i ) {
			String chr = Integer.toString(i);
			genomeLength += maizeGenome.get(chr).length();
		}
		
		HashMap<String, ArrayList<GeneSimple>> geneHashMap = ifIntron.getGeneHashMap();
		HashMap<String, HashSet<Integer>> allGeneticsBps = new HashMap<String, HashSet<Integer>>();
		int allGeneticLength=0;
		for( int i=1; i<=10; i++ ) {
			String chr = Integer.toString(i);
			allGeneticsBps.put(chr, new HashSet<Integer>());
			for( GeneSimple geneSimple : geneHashMap.get(chr) ) {
				for( int position=geneSimple.getStart(); position<=geneSimple.getEnd(); ++position ) {
					allGeneticsBps.get(chr).add(position);
				}
			}
			allGeneticLength += allGeneticsBps.get(chr).size();
		}
		
		
		
		HashSet<String> mCG = new HashSet<String>();
		mCG.add("CGA");
		mCG.add("CGC");
		mCG.add("CGG");
		mCG.add("CGT");
		mCG.add("CGN");
		
		HashSet<String> CHG = new HashSet<String>();
		CHG.add("CAG");
		CHG.add("CCG");
		CHG.add("CTG");
		
		HashSet<String> CHH = new HashSet<String>();
		CHH.add("CAA");
		CHH.add("CCA");
		CHH.add("CTA");
		CHH.add("CAC");
		CHH.add("CCC");
		CHH.add("CTC");
		CHH.add("CAT");
		CHH.add("CCT");
		CHH.add("CTT");
		
		try {
			HashMap<String, HashSet<Integer>> nonIntronUtrCns = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> nonOpTfbsEqtlIntronUtrCns = new HashMap<String, HashSet<Integer>> ();
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				nonIntronUtrCns.put(chr, new HashSet<Integer>());
				nonOpTfbsEqtlIntronUtrCns.put(chr, new HashSet<Integer>());
			}
			readWigFile(nonIntronUtrCnsWigFile, nonIntronUtrCns);
			readWigFile(nonOpTfbsEqtlIntronUtrCnsWigFile, nonOpTfbsEqtlIntronUtrCns);
			
			try {
				int otherNumberOfCGSites = 0;
				int otherNumberOfCHGSites = 0;
				int otherNumberOfCHHSites = 0;
				
				int genicNumberOfCGSites = 0;
				int genicNumberOfCHGSites = 0;
				int genicNumberOfCHHSites = 0;
				
				int nonIntronUtrCnsNumberOfCGSites = 0;
				int nonIntronUtrCnsNumberOfCHGSites = 0;
				int nonIntronUtrCnsNumberOfCHHSites = 0;
				
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCGSites = 0;
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCHGSites = 0;
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCHHSites = 0;
				
				
				int otherNumberOfCGMethylations = 0;
				int otherNumberOfCHGMethylations = 0;
				int otherNumberOfCHHMethylations = 0;
				
				int genicNumberOfCGMethylations = 0;
				int genicNumberOfCHGMethylations = 0;
				int genicNumberOfCHHMethylations = 0;
				
				int nonIntronUtrCnsNumberOfCGMethylations = 0;
				int nonIntronUtrCnsNumberOfCHGMethylations = 0;
				int nonIntronUtrCnsNumberOfCHHMethylations = 0;
				
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCGMethylations = 0;
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCHGMethylations = 0;
				int nonOpTfbsEqtlIntronUtrCnsNumberOfCHHMethylations = 0;
				
				int nonIntronUtrCnsLength = 0;
				int nonOpTfbsEqtlIntronUtrCnsLength = 0;
				for( int i=1; i<=10; ++i ) {
					String chr = Integer.toString(i);
					nonIntronUtrCnsLength += nonIntronUtrCns.get(chr).size();
					nonOpTfbsEqtlIntronUtrCnsLength += nonOpTfbsEqtlIntronUtrCns.get(chr).size();
				}
				File file = new File( "/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_allc_methylome_B73_leaf.tsv" );
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = null;
				while ((tempString = reader.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+");
					String chr = arrOfStr[0];
					if( nonIntronUtrCns.containsKey(chr) ) {
						String seq = arrOfStr[3];
						
						int position = Integer.parseInt( arrOfStr[1] );
						boolean methylated = false;
						if( arrOfStr[6].charAt(0) == '1' ) {
							methylated = true;
						}
						if( mCG.contains(seq) ) {
							if( ifIntron.getElement(chr, position) == 0 ) {
								if( nonOpTfbsEqtlIntronUtrCns.get(chr).contains(position) ) {
									++nonOpTfbsEqtlIntronUtrCnsNumberOfCGSites;
									if( methylated ) {
										++nonOpTfbsEqtlIntronUtrCnsNumberOfCGMethylations;
									}
								}else if( nonIntronUtrCns.get(chr).contains(position) ) { // if it is not unknown, then it is known
									++nonIntronUtrCnsNumberOfCGSites;
									if( methylated ) {
										++nonIntronUtrCnsNumberOfCGMethylations;
									}
								} else {
									++otherNumberOfCGSites;
									if( methylated ) {
										++otherNumberOfCGMethylations;
									}
								}
							}else {
								++genicNumberOfCGSites;
								if( methylated ) {
									++genicNumberOfCGMethylations;
								}
							}
						}else if( CHG.contains(seq) ) {
							if( ifIntron.getElement(chr, position) == 0 ) {
								if( nonOpTfbsEqtlIntronUtrCns.get(chr).contains(position) ) {
									++nonOpTfbsEqtlIntronUtrCnsNumberOfCHGSites;
									if( methylated ) {
										++nonOpTfbsEqtlIntronUtrCnsNumberOfCHGMethylations;
									}
								}else if( nonIntronUtrCns.get(chr).contains(position) ) {
									++nonIntronUtrCnsNumberOfCHGSites;
									if( methylated ) {
										++nonIntronUtrCnsNumberOfCHGMethylations;
									}
								}else {
									++otherNumberOfCHGSites;
									if( methylated ) {
										++otherNumberOfCHGMethylations;
									}
								}
							}else {
								++genicNumberOfCHGSites;
								if( methylated ) {
									++genicNumberOfCHGMethylations;
								}
							}
						}else if( CHH.contains(seq) ) {
							if( ifIntron.getElement(chr, position) == 0 ) {
								if( nonOpTfbsEqtlIntronUtrCns.get(chr).contains(position) ) {
									++nonOpTfbsEqtlIntronUtrCnsNumberOfCHHSites;
									if( methylated ) {
										++nonOpTfbsEqtlIntronUtrCnsNumberOfCHHMethylations;
									}
								}else if( nonIntronUtrCns.get(chr).contains(position) ) {
									++nonIntronUtrCnsNumberOfCHHSites;
									if( methylated ) {
										++nonIntronUtrCnsNumberOfCHHMethylations;
									}
								}else {
									++otherNumberOfCHHSites;
									if( methylated ) {
										++otherNumberOfCHHMethylations;
									}
								}
							}else {
								++genicNumberOfCHHSites;
								if( methylated ) {
									++genicNumberOfCHHMethylations;
								}
							}
						}
					}
				}
				reader.close();
				
				System.out.println("otherNumberOfCGSites=" + otherNumberOfCGSites + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCGSites=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCGSites + "\nnonIntronUtrCnsNumberOfCGSites=" + nonIntronUtrCnsNumberOfCGSites + "\ngenicNumberOfCGSites=" + genicNumberOfCGSites);
				System.out.println("otherNumberOfCHGSites=" + otherNumberOfCHGSites + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCHGSites=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCHGSites + "\nnonIntronUtrCnsNumberOfCHGSites=" + nonIntronUtrCnsNumberOfCHGSites + "\ngenicNumberOfCHGSites=" + genicNumberOfCHGSites);
				System.out.println("otherNumberOfCHHSites=" + otherNumberOfCHHSites + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCHHSites=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCHHSites + "\nnonIntronUtrCnsNumberOfCHHSites=" + nonIntronUtrCnsNumberOfCHHSites + "\ngenicNumberOfCHHSites=" + genicNumberOfCHHSites);
				
				System.out.println("otherNumberOfCGMethylations=" + otherNumberOfCGMethylations + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCGMethylations=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCGMethylations + "\nnonIntronUtrCnsNumberOfCGMethylations=" + nonIntronUtrCnsNumberOfCGMethylations + "\ngenicNumberOfCGMethylations=" + genicNumberOfCGMethylations);
				System.out.println("otherNumberOfCHGMethylations=" + otherNumberOfCHGMethylations + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCHGMethylations=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCHGMethylations + "\nnonIntronUtrCnsNumberOfCHGMethylations=" + nonIntronUtrCnsNumberOfCHGMethylations + "\ngenicNumberOfCHGMethylations=" + genicNumberOfCHGMethylations);
				System.out.println("otherNumberOfCHHMethylations=" + otherNumberOfCHHMethylations + 
						"\nnonOpTfbsEqtlIntronUtrCnsNumberOfCHHMethylations=" + nonOpTfbsEqtlIntronUtrCnsNumberOfCHHMethylations + "\nnonIntronUtrCnsNumberOfCHHMethylations=" + nonIntronUtrCnsNumberOfCHHMethylations + "\ngenicNumberOfCHHMethylations=" + genicNumberOfCHHMethylations);
				nonIntronUtrCnsLength = nonIntronUtrCnsLength - nonOpTfbsEqtlIntronUtrCnsLength;  // if it is not unknown, then it is known
				System.out.println("nonIntronUtrCnsLength=" +  nonIntronUtrCnsLength);
				System.out.println("nonOpTfbsEqtlIntronUtrCnsLength=" + nonOpTfbsEqtlIntronUtrCnsLength);
				System.out.println("allGeneticLength=" + allGeneticLength);
				int otherSeqLength = genomeLength - nonIntronUtrCnsLength - nonOpTfbsEqtlIntronUtrCnsLength - allGeneticLength;
				System.out.println("otherSeqLength=" + otherSeqLength);
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public OverLapCnsWithMethalation(String wigFile) {
		
		int genomeLength = 0;
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		for( int i=1; i<=10; ++i ) {
			String chr = Integer.toString(i);
			genomeLength += maizeGenome.get(chr).length();
		}
		
		HashSet<String> mCG = new HashSet<String>();
		mCG.add("CGA");
		mCG.add("CGC");
		mCG.add("CGG");
		mCG.add("CGT");
		mCG.add("CGN");
		
		HashSet<String> CHG = new HashSet<String>();
		CHG.add("CAG");
		CHG.add("CCG");
		CHG.add("CTG");
		
		HashSet<String> CHH = new HashSet<String>();
		CHH.add("CAA");
		CHH.add("CCA");
		CHH.add("CTA");
		CHH.add("CAC");
		CHH.add("CCC");
		CHH.add("CTC");
		CHH.add("CAT");
		CHH.add("CCT");
		CHH.add("CTT");
		
		try {
			HashMap<String, HashSet<Integer>> cns = new HashMap<String, HashSet<Integer>> ();
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				cns.put(chr, new HashSet<Integer>());
			}
			readWigFile(wigFile, cns);
			
			try {
				int otherNumberOfCGSites = 0;
				int otherNumberOfCHGSites = 0;
				int otherNumberOfCHHSites = 0;
				
				int cnsNumberOfCGSites = 0;
				int cnsNumberOfCHGSites = 0;
				int cnsNumberOfCHHSites = 0;
				
				int otherNumberOfCGMethylations = 0;
				int otherNumberOfCHGMethylations = 0;
				int otherNumberOfCHHMethylations = 0;
				
				int cnsNumberOfCGMethylations = 0;
				int cnsNumberOfCHGMethylations = 0;
				int cnsNumberOfCHHMethylations = 0;
				
				int cnsLength = 0;
				for( int i=1; i<=10; ++i ) {
					String chr = Integer.toString(i);
					cnsLength += cns.get(chr).size();
				}
				File file = new File( "/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_allc_methylome_B73_leaf.tsv" );
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = null;
				while ((tempString = reader.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+");
					String chr = arrOfStr[0];
					if( cns.containsKey(chr) ) {
						String seq = arrOfStr[3];
						
						int position = Integer.parseInt( arrOfStr[1] );
						boolean methylated = false;
						if( arrOfStr[6].charAt(0) == '1' ) {
							methylated = true;
						}
						if( mCG.contains(seq) ) {
							if( cns.get(chr).contains(position) ) {
								++cnsNumberOfCGSites;
								if( methylated ) {
									++cnsNumberOfCGMethylations;
								}
							}else {
								++otherNumberOfCGSites;
								if( methylated ) {
									++otherNumberOfCGMethylations;
								}
							}
						}else if( CHG.contains(seq) ) {
							if( cns.get(chr).contains(position) ) {
								++cnsNumberOfCHGSites;
								if( methylated ) {
									++cnsNumberOfCHGMethylations;
								}
							}else {
								++otherNumberOfCHGSites;
								if( methylated ) {
									++otherNumberOfCHGMethylations;
								}
							}
						
						}else if( CHH.contains(seq) ) {
							if( cns.get(chr).contains(position) ) {
								++cnsNumberOfCHHSites;
								if( methylated ) {
									++cnsNumberOfCHHMethylations;
								}
							}else {
								++otherNumberOfCHHSites;
								if( methylated ) {
									++otherNumberOfCHHMethylations;
								}
							}
						
						}
					}
				}
				reader.close();
				
				System.out.println("otherNumberOfCGSites=" + otherNumberOfCGSites + 
						"\ncnsNumberOfCGSites=" + cnsNumberOfCGSites );
				System.out.println("otherNumberOfCHGSites=" + otherNumberOfCHGSites + 
						"\ncnsNumberOfCHGSites=" + cnsNumberOfCHGSites );
				System.out.println("otherNumberOfCHHSites=" + otherNumberOfCHHSites + 
						"\ncnsNumberOfCHHSites=" + cnsNumberOfCHHSites);
				
				System.out.println("otherNumberOfCGMethylations=" + otherNumberOfCGMethylations + 
						"\ncnsNumberOfCGMethylations=" + cnsNumberOfCGMethylations );
				System.out.println("otherNumberOfCHGMethylations=" + otherNumberOfCHGMethylations + 
						"\ncnsNumberOfCHGMethylations=" + cnsNumberOfCHGMethylations);
				System.out.println("otherNumberOfCHHMethylations=" + otherNumberOfCHHMethylations + 
						"\ncnsNumberOfCHHMethylations=" + cnsNumberOfCHHMethylations );
				System.out.println("cnsLength=" +  cnsLength);
				int otherSeqLength = genomeLength - cnsLength;
				System.out.println("otherSeqLength=" + otherSeqLength);
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void readWigFile(String wigFile, HashMap<String, HashSet<Integer>> wigRecords) {
		try {
			File file = new File( wigFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String chr = "";
			Pattern p =Pattern.compile("chrom=(\\S+)");
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) == 't' ) {
					
				}else if( tempString.charAt(0) == 'v' ) {
					Matcher m = p.matcher(tempString);
					if( m.find() ) {
						chr = m.group(1);
					}
				}else {
					String[] arrOfStr = tempString.split("\\s+");
					if( wigRecords.containsKey(chr) ) {
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							wigRecords.get(chr).add(Integer.parseInt(arrOfStr[0]));
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
