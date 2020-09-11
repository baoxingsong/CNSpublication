package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.LeadingSnps;

public class OverLapCnsWitLeadingEqtl {
	public static void main(String[] args) {
		// only use those CNS never overlap genetic region
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		
		new OverLapCnsWitLeadingEqtl("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", maizeGenome);
		new OverLapCnsWitLeadingEqtl("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBSOrEnhancer.wig", maizeGenome);
		new OverLapCnsWitLeadingEqtl("/media/bs674/1_8t/AndCns/CNSprofiles/loop.wig", maizeGenome);
	}
	
	public OverLapCnsWitLeadingEqtl(String wigFile, HashMap<String, String> maizeGenome) {
		System.out.println("wigFile:" + wigFile);
		try {
			LeadingSnps leadingSnps = new LeadingSnps();
			leadingSnps.readSummary();
			HashMap<String, HashSet<Integer>> leadingEqtl = leadingSnps.getLeadingSnps();
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			
			readWigFile(wigFile, cnsRecords);
			System.out.println("bed files reading done");
			
			HashMap<Integer, Integer> classes = new HashMap<Integer, Integer>();
    		for( String chr : cnsRecords.keySet() ) { // so here only chromosomes 1-10 were used
    			for( int i=0; i<maizeGenome.get(chr).length(); ++i ) {
    				int thisValue = 0;
					if( cnsRecords.get(chr).contains(i+1)  ) {
    					thisValue =  (thisValue + 1);
    				}
    				
    				if( leadingEqtl.get(chr).contains(i+1)  ) {
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
