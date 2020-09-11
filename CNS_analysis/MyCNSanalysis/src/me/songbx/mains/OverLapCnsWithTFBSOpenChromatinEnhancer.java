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
import me.songbx.model.Cds;
import me.songbx.model.GeneSimple;
import me.songbx.model.Transcript;
import me.songbx.service.IfIntron;

public class OverLapCnsWithTFBSOpenChromatinEnhancer {
	public static void main(String[] args) {
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		HashMap<String, ArrayList<Transcript>> transcriptHashMap = ifIntron.getTranscriptHashMap();
		HashMap<String, ArrayList<GeneSimple>> geneHashMap = ifIntron.getGeneHashMap();
		HashMap<String, HashSet<Integer>> allGeneticsBps = new HashMap<String, HashSet<Integer>>();
		HashMap<String, HashSet<Integer>> allCdsBps = new HashMap<String, HashSet<Integer>>();
		for( int i=1; i<=10; i++ ) {
			String chr = Integer.toString(i);
			allGeneticsBps.put(chr, new HashSet<Integer>());
			for( GeneSimple geneSimple : geneHashMap.get(chr) ) {
				for( int position=geneSimple.getStart(); position<=geneSimple.getEnd(); ++position ) {
					allGeneticsBps.get(chr).add(position);
				}
			}
			allCdsBps.put(chr, new HashSet<Integer>());
			for( Transcript transcript : transcriptHashMap.get(chr) ) {
				for( Cds cds :  transcript.getCdsHashSet()) {
					for( int position=cds.getStart(); position<=cds.getEnd(); ++position ) {
						allCdsBps.get(chr).add(position);
					}
				}
			}
		}
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.wig", allCdsBps);
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.wig", allCdsBps);
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_depth.wig", allCdsBps);
		
		// if you only want to check the intergenic region
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.wig", allGeneticsBps);
		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.wig", allGeneticsBps);
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_depth.wig", allGeneticsBps);
//		new OverLapCnsWithTFBSOpenChromatinEnhancer("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", allGeneticsBps);
    }
	
	public OverLapCnsWithTFBSOpenChromatinEnhancer(String wigFile, HashMap<String, HashSet<Integer>> excludeBps) {
		try {
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
			HashMap<String, HashSet<Integer>> tFBSBedRecords = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> openChromatinRecords = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> enhancers = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				tFBSBedRecords.put(chr, new HashSet<Integer>());
				openChromatinRecords.put(chr, new HashSet<Integer>());
				enhancers.put(chr, new HashSet<Integer>());
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			
			String tfbs_bed_file = "/media/bs674/2t/openChromotainAndTFBS/all_reproducible_peaks_summits_merged.bed/all_reproducible_peaks_summits_merged.bed";
			readBedFile(tfbs_bed_file, tFBSBedRecords);
			
			ArrayList<String> enhancerFiles = new ArrayList<String>();
			enhancerFiles.add("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed");
			enhancerFiles.add("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed");
			
			
			ArrayList<String> openChromatinFiles = new ArrayList<String>();
			openChromatinFiles.add("/media/bs674/2t/openChromotainAndTFBS/ATAC_seq_data/BEDfiles/Zm_anotation_rep1.intersect.rep2.acr.bed");
			openChromatinFiles.add("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed");
			openChromatinFiles.add("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398047_ATAC_B73_ear.filtered_ACR.bed");
			for( String openChromatinFile : openChromatinFiles ) {
				readBedFile(openChromatinFile, openChromatinRecords);
			}
			
			for( String enhancerFile : enhancerFiles ) {
				readBedFile(enhancerFile, enhancers);
			}
			readWigFile(wigFile, cnsRecords);
			System.out.println("bed files reading done");
			
			HashMap<Short, Integer> classes = new HashMap<Short, Integer>();
    		for( String chr : tFBSBedRecords.keySet() ) { // so here only chromosomes 1-10 were used
    			for( int i=0; i<maizeGenome.get(chr).length(); ++i ) {
    				short thisValue = 0;
    				if( excludeBps.get(chr).contains(i+1) ) {
    					thisValue=-1;
    				}else {
						if( cnsRecords.get(chr).contains(i+1)  ) {
	    					thisValue = (short) (thisValue + 1);
	    				}
	    				if( enhancers.containsKey(chr) && enhancers.get(chr).contains(i+1) ) {
	    					thisValue = (short) (thisValue + 2);
	    				}
	    				if( tFBSBedRecords.get(chr).contains(i+1)  ) {
	    					thisValue = (short) (thisValue + 4);
	    				}
	    				if( openChromatinRecords.get(chr).contains(i+1) ) {
	    					thisValue = (short) (thisValue + 8);
	    				}
    				}
    					
    				if( classes.containsKey(thisValue) ) {
    					classes.put(thisValue, classes.get(thisValue)+1);
    				}else {
    					classes.put(thisValue, 1);
    				}
    			}
    		}
    		System.out.println(wigFile + " done");
    		for( short key : classes.keySet() ) {
    			System.out.println("n"+key + " = " + classes.get(key));
    		}
    		
    		try {
    			int totalRecords = 0;
    			int touchedRecords = 0;
    			File file = new File( tfbs_bed_file );
    			BufferedReader reader = new BufferedReader(new FileReader(file));
    			String tempString = null;
    			while ((tempString = reader.readLine()) != null) {
    				String[] arrOfStr = tempString.split("\\s+");
    				if( cnsRecords.containsKey(arrOfStr[0]) ) {
    					++totalRecords;
    					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
    						if( cnsRecords.get(arrOfStr[0]).contains(i)  ) {
    							++touchedRecords;
    							break;
    						}
    					}
    				}
    			}
    			reader.close();
    			System.out.println("total TBFS records:" + totalRecords + " total touched by cns:" + touchedRecords);
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
    		
    		try {
    			int totalRecords = 0;
    			int touchedRecords = 0;
	    		for( String openChromatinFile : openChromatinFiles ) {
	    			int thisTotalRecords = 0;
	    			int thisTouchedRecords = 0;
	    			File file = new File( openChromatinFile );
	    			BufferedReader reader = new BufferedReader(new FileReader(file));
	    			String tempString = null;
	    			while ((tempString = reader.readLine()) != null) {
	    				String[] arrOfStr = tempString.split("\\s+");
	    				if( cnsRecords.containsKey(arrOfStr[0]) ) {
	    					++totalRecords;
	    					++thisTotalRecords;
	    					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
	    						if( cnsRecords.get(arrOfStr[0]).contains(i)  ) {
	    							++touchedRecords;
	    							++thisTouchedRecords;
	    							break;
	    						}
	    					}
	    				}
	    			}
	    			reader.close();
	    			System.out.println(openChromatinFile + " open-Chromatin records:" + thisTotalRecords + " touched by cns:" + thisTouchedRecords);
    			}
    			System.out.println("total open-Chromatin records:" + totalRecords + " total touched by cns:" + touchedRecords);
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
    		
    		
    		try {
    			int totalRecords = 0;
    			int touchedRecords = 0;
    			ArrayList<String> openChromatinFiles2 = new ArrayList<String>();
    			openChromatinFiles2.add("/media/bs674/1_8t/AndCns/openchromatin/AP.bfthresh1.1.MNaseHS.Ranges.v4.bed");
    			openChromatinFiles2.add("/media/bs674/1_8t/AndCns/openchromatin/RP.bfthresh1.1.MNaseHS.Ranges.v4.bed");
	    		for( String openChromatinFile : openChromatinFiles2 ) {
	    			int thisTotalRecords = 0;
	    			int thisTouchedRecords = 0;
	    			File file = new File( openChromatinFile );
	    			BufferedReader reader = new BufferedReader(new FileReader(file));
	    			String tempString = null;
	    			while ((tempString = reader.readLine()) != null) {
	    				String[] arrOfStr = tempString.split("\\s+");
	    				if( cnsRecords.containsKey(arrOfStr[0]) ) {
	    					++totalRecords;
	    					++thisTotalRecords;
	    					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
	    						if( cnsRecords.get(arrOfStr[0]).contains(i)  ) {
	    							++touchedRecords;
	    							++thisTouchedRecords;
	    							break;
	    						}
	    					}
	    				}
	    			}
	    			reader.close();
	    			System.out.println(openChromatinFile + " open-Chromatin records:" + thisTotalRecords + " touched by cns:" + thisTouchedRecords);
    			}
    			System.out.println("total open-Chromatin records:" + totalRecords + " total touched by cns:" + touchedRecords);
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
    		
    		try {
    			int totalRecords = 0;
    			int touchedRecords = 0;
    			for( String enhancerFile : enhancerFiles ) {
    				int thisTotalRecords = 0;
	    			int thisTouchedRecords = 0;
    				File file = new File( enhancerFile );
	    			BufferedReader reader = new BufferedReader(new FileReader(file));
	    			String tempString = null;
	    			while ((tempString = reader.readLine()) != null) {
	    				String[] arrOfStr = tempString.split("\\s+");
	    				if( cnsRecords.containsKey(arrOfStr[0]) ) {
	    					++totalRecords;
	    					++thisTotalRecords;
	    					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
	    						if( cnsRecords.get(arrOfStr[0]).contains(i)  ) {
	    							++touchedRecords;
	    							++thisTouchedRecords;
	    							break;
	    						}
	    					}
	    				}
	    			}
	    			reader.close();
	    			System.out.println(enhancerFile + " open-Chromatin records:" + thisTotalRecords + " touched by cns:" + thisTouchedRecords );
    			}
    			System.out.println("total enhancer records:" + totalRecords + " total touched by cns:" + touchedRecords);
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
    		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public void readBedFile(String bedFile, HashMap<String, HashSet<Integer>> bedRecords) {
		try {
			File file = new File( bedFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				if( bedRecords.containsKey(arrOfStr[0]) ) {
					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
						bedRecords.get(arrOfStr[0]).add(i);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
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
