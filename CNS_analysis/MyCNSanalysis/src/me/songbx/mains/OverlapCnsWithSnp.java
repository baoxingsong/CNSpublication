package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import me.songbx.impl.ChromoSomeReadImpl;
import htsjdk.tribble.readers.TabixReader;

public class OverlapCnsWithSnp {
	public static void main(String[] args) {
		try {
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/Zea_mays.AGPv4.dna.toplevel_callable.fa");
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
			HashMap<String, ArrayList<Integer>> callableCoverage = new HashMap<String, ArrayList<Integer>>();
			for( String chr : maizeGenome.keySet() ) {
				callableCoverage.put( chr, new ArrayList<Integer>() );
				for (int i = 0; i < maizeGenome.get(chr).length(); ++i) {
					callableCoverage.get(chr).add(0);
				}
			}
			
        	File file = new File("/media/bs674/and-CNS/maize282Genotyping/callablePipeline/callable_list");
        	int coutOfCoverageFiles = 0;
        	BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            while ((tempString = reader.readLine()) != null) {
            	tempString = tempString.trim();
            	if( tempString.length() > 3 ) {
            		++coutOfCoverageFiles;
            	}
            }
            reader.close();
            //int miniCoverate = coutOfCoverageFiles/2;
            int miniCoverate = 50;
            
        	reader = new BufferedReader(new FileReader(file));
            while ((tempString = reader.readLine()) != null) {
            	tempString = tempString.trim();
            	if( tempString.length() > 3 ) {
	            	System.out.println(tempString);
	        		TabixReader tr = new TabixReader(tempString);
	        		String s;
					while ((s = tr.readLine() ) != null) {
						String[] arrOfStr = s.split("\\s+"); 
						if( arrOfStr[3].compareTo("CALLABLE")==0 ) {
							for( int pi=Integer.parseInt(arrOfStr[1]); pi< Integer.parseInt(arrOfStr[2]); ++pi ) {
								if( callableCoverage.get(arrOfStr[0]).get(pi) < miniCoverate) {
									callableCoverage.get(arrOfStr[0]).set(pi, callableCoverage.get(arrOfStr[0]).get(pi) + 1);
								}
							}
						}
						arrOfStr = null;
					}
					tr.close();
					System.gc();
            	}
            }
            reader.close();
            
            for( String chr : maizeGenome.keySet() ) {
            	outPut.write(">" + chr + "\n");
            	for (int i = 0; i < maizeGenome.get(chr).length(); i++) {
            		if ( callableCoverage.get(chr).get(i)>= miniCoverate) {
                		outPut.write(maizeGenome.get(chr).charAt(i));            			
            		}else {
            			outPut.write("b"); // code no-callable as b     		
            		}
            		if ( 59 == i%60  ) {
            			outPut.write("\n");
            		}
            	}
            	outPut.write("\n");
            }
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
