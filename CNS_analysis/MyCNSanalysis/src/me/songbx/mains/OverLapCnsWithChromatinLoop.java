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

public class OverLapCnsWithChromatinLoop {
	public static void main(String[] args) {
		// only use those CNS never overlap genetic region
		
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1twdblack/mingqiuAndHai_CNS/maizev4/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		
		new OverLapCnsWithChromatinLoop("/media/bs674/wd/CNSprofiles/nonIntronUtrCns.wig", maizeGenome);
		new OverLapCnsWithChromatinLoop("/media/bs674/wd/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.wig", maizeGenome);
		
	}
	
	public OverLapCnsWithChromatinLoop(String wigFile, HashMap<String, String> maizeGenome) {
		System.out.println("wigFile:" + wigFile);
		try {
			HashMap<String, HashSet<Integer>> chromationLoopRecords = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				chromationLoopRecords.put(chr, new HashSet<Integer>());
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			
			ArrayList<String> chromation_loops = new ArrayList<String>();
			chromation_loops.add("/media/bs674/wd/maizeHicAndMethylation/GSE120304_HiC_B73_leaf_merged_loops.txt");
			chromation_loops.add("/media/bs674/wd/maizeHicAndMethylation/RNAPII.loop");
			chromation_loops.add("/media/bs674/wd/maizeHicAndMethylation/h3k4me3.loop");
			
			//chromation_loops.add("/media/bs674/wd/maizeHicAndMethylation/Corteva/stable1");
			//chromation_loops.add("/media/bs674/wd/maizeHicAndMethylation/Corteva/stable2");
			for(String chromation_loop : chromation_loops) {
				readBedFile(chromation_loop, chromationLoopRecords);
			}
			readWigFile(wigFile, cnsRecords);
			System.out.println("bed files reading done");
			
			HashMap<Integer, Integer> classes = new HashMap<Integer, Integer>();
    		for( String chr : chromationLoopRecords.keySet() ) { // so here only chromosomes 1-10 were used
    			for( int i=0; i<maizeGenome.get(chr).length(); ++i ) {
    				int thisValue = 0;
					if( cnsRecords.get(chr).contains(i+1)  ) {
    					thisValue =  (thisValue + 1);
    				}
    				
    				if( chromationLoopRecords.get(chr).contains(i+1)  ) {
    					thisValue = (thisValue + 2);
    				}
    				
    				if( classes.containsKey(thisValue) ) {
    					classes.put(thisValue, classes.get(thisValue)+1);
    				}else {
    					classes.put(thisValue, 1);
    				}
    			}
    		}
 
    		for( int key : classes.keySet() ) {
    			System.out.println("n"+key + " = " + classes.get(key));
    		}
    		int totalLength = classes.get(0) +  classes.get(1) + classes.get(2) + classes.get(3);
    		double enrichment = ((double)classes.get(3)/(double)classes.get(2) )/((double)(classes.get(1)+classes.get(3))/(double)totalLength) ;
    		System.out.println("enrichment:" + enrichment);
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
				//if ( arrOfStr[0].charAt(0) != arrOfStr[3].charAt(0) ) {
					if( bedRecords.containsKey(arrOfStr[0]) ) {
						for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
							bedRecords.get(arrOfStr[0]).add(i);
						}
					}
					if( bedRecords.containsKey(arrOfStr[3]) ) {
						for( int i=Integer.parseInt(arrOfStr[4])+1; i<=Integer.parseInt(arrOfStr[5]); ++i ) {
							bedRecords.get(arrOfStr[3]).add(i);
						}
					}
				//}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void readBedFile2(String bedFile, HashMap<String, HashSet<Integer>> bedRecords) {
		try {
			File file = new File( bedFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				//if ( arrOfStr[0].charAt(0) != arrOfStr[3].charAt(0) ) {
					if( bedRecords.containsKey(arrOfStr[0]) ) {
						for( int i=Integer.parseInt(arrOfStr[1])-10000; i<=Integer.parseInt(arrOfStr[1]); ++i ) {
							bedRecords.get(arrOfStr[0]).add(i);
						}
						for( int i=Integer.parseInt(arrOfStr[2])+1; i<=Integer.parseInt(arrOfStr[2])+10000; ++i ) {
							bedRecords.get(arrOfStr[0]).add(i);
						}
					}
					if( bedRecords.containsKey(arrOfStr[3]) ) {
						for( int i=Integer.parseInt(arrOfStr[4])-10000; i<=Integer.parseInt(arrOfStr[4]); ++i ) {
							bedRecords.get(arrOfStr[3]).add(i);
						}
						for( int i=Integer.parseInt(arrOfStr[5])+1; i<=Integer.parseInt(arrOfStr[5])+10000; ++i ) {
							bedRecords.get(arrOfStr[3]).add(i);
						}
					}
				//}
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
