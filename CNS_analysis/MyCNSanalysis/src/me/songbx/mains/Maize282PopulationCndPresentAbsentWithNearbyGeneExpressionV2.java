package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.model.BedRecord;
import htsjdk.tribble.readers.TabixReader;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;

import java.nio.file.Path;
import java.nio.file.Paths;


public class Maize282PopulationCndPresentAbsentWithNearbyGeneExpressionV2 {
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSBasedGwas/genotype");
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			HashMap<String, WigFileReader> wigFileReaders = new HashMap<String, WigFileReader>();
			HashMap<String, WigFileReader> wigFileReaders2 = new HashMap<String, WigFileReader>();
			
			outPut.print("GenetypeID\tlength");
			HashMap<String, String> seqNameToId = new HashMap<String, String>();
	        try {
	        	File fastaFile = new File("/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/cns_seq.fa");
	        	BufferedReader fastaReader = new BufferedReader(new FileReader(fastaFile));
	        	String tempString = null;
	        	Pattern pattern = Pattern.compile("^>(\\S+)\\s+(\\S+)");
	        	while ((tempString = fastaReader.readLine()) != null) {
	        		Matcher match = pattern.matcher(tempString);
	            	if( match.find() ) {
	            		seqNameToId.put(match.group(2), match.group(1));
	            	}
	        	}
	        	fastaReader.close();
	        	
	        	File file = new File(args[0]);
	        	BufferedReader reader = new BufferedReader(new FileReader(file));
	            
	            while ((tempString = reader.readLine()) != null) {
	            	tempString = tempString.trim();
	        		String baseName = tempString;
	        		baseName = baseName.replaceAll(".*\\/", "");
	        		baseName = baseName.replaceAll("_callable_status.bed.gz", "");
	        		outPut.print("\t" + baseName);
	        		
	        		String bwFilePosition0 =   "/media/bs674/panAndAssemblyfi/maize282Genotyping/coverageBwFiles/" + baseName  + ".bw";
            		Path bwFile0 = Paths.get(bwFilePosition0);
            		WigFileReader wig0 = WigFileReader.autodetect(bwFile0);
            		wigFileReaders.put(baseName, wig0);
            		
            		String bwFilePosition =   "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/mapreadsToCns/" + baseName  + ".bw";
            		Path bwFile = Paths.get(bwFilePosition);
            		WigFileReader wig = WigFileReader.autodetect(bwFile);
            		wigFileReaders2.put(baseName, wig);
            		
					TabixReader tr = new TabixReader(tempString);
					tabixReaders.put(baseName, tr);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        outPut.println();
	        
	        System.out.println("begin msa file reading");
	        HashSet<String> titleLines = new HashSet<String>();
	        File file0 = new File("/media/bs674/1_8t/AndCns/msa");
	        BufferedReader reader0 = new BufferedReader(new FileReader(file0));
            String tempString0 = null;
            while ((tempString0 = reader0.readLine()) != null) {
            	if( tempString0.startsWith(">reference") ) {
            		titleLines.add(tempString0);
            	}
            }
            reader0.close();
            System.out.println("msa file reading done");
            
            Pattern p = Pattern.compile("^>reference:\\s+(\\S+):(\\d+)-(\\d+)$");
            for (String tempString1 : titleLines) {
            	Matcher m=p.matcher(tempString1);
            	if( m.find() ) {
            		int start = Integer.parseInt(m.group(2));
    	            int end = Integer.parseInt(m.group(3));
    	            String chr = m.group(1);
    	            HashSet<Integer> nonMaskingPositions = new HashSet<Integer>();
    	            int totalLength = 0;
    	            for( int i=start-1; i<end; ++i ) {
    	            	if( maizeGenome.get(chr).charAt(i) != 'n' ) {
            				totalLength = totalLength + 1;
            				nonMaskingPositions.add(i);
    	            	}
    	            }
    	            if(totalLength > 30) {
    		            outPut.print(chr+":"+start+"-"+end+"\t" + totalLength);
    		            String seq_name = chr+":"+start+"-"+end;
    		            
    		            //read bed files begin
    		    		File file = new File(args[0]);
    		            try {
    		            	BufferedReader reader = new BufferedReader(new FileReader(file));
    		                String tempString = null;
    		                while ((tempString = reader.readLine()) != null) {
    		                	tempString = tempString.trim();
    		                	try {
    		                		ArrayList<BedRecord> bedRecords = new ArrayList<BedRecord>();
    		                		String baseName = tempString;
    		                		baseName = baseName.replaceAll(".*\\/", "");
    		                		baseName = baseName.replaceAll("_callable_status.bed.gz", "");
    		                		WigFileReader wig = wigFileReaders.get(baseName);
    		                		WigFileReader wig2 = wigFileReaders2.get(baseName);
    		                		TabixReader tr = tabixReaders.get(baseName);
    	
    		                		String s;
    		    					TabixReader.Iterator iter = tr.query(chr+":"+start+"-"+end); // get the iterator
    		    					while ((s = iter.next()) != null) {
    		    						String[] arrOfStr = s.split("\\s+"); 
    		    						if( arrOfStr[3].compareTo("CALLABLE")==0 ) {
    		    							BedRecord bedRecord = new BedRecord(Integer.parseInt(arrOfStr[1])+1, Integer.parseInt(arrOfStr[2])); //bed has special coordinate
    		    							bedRecords.add(bedRecord);
    		    						}
    		    					}
    		    					int totalGoodLength = 0;
    		    		            for( int i : nonMaskingPositions ) {
    		    		            	//if( maizeGenome.get(chr).charAt(i) != 'n' && ifCallAble(i+1, bedRecords, wig, chr, i+2-start, wig2, seqNameToId.get(seq_name)) ) {
    		    		            	if( maizeGenome.get(chr).charAt(i) != 'n' && ifCallAble(i+1, bedRecords, wig, chr) ) {
    		    		            		totalGoodLength = totalGoodLength + 1;
    		    		            	}
    		    		            }
    		    		            outPut.print("\t" + totalGoodLength);
    	
    		                	} catch (IOException e) {
    		    					e.printStackTrace();
    		    				}
    		                }
    		                reader.close();
    		                outPut.println();
    		            } catch (IOException e) {
    		                e.printStackTrace();
    		            }//read bed files end
    	            }
            	}
			}
            outPut.close();
			for( String baseName : tabixReaders.keySet()  ) {
				tabixReaders.get(baseName).close();
				wigFileReaders.get(baseName).close();
				wigFileReaders2.get(baseName).close();
			}
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
	
	static boolean ifCallAble( int position, ArrayList<BedRecord> bedRecords, WigFileReader wig, String chr, int position2, WigFileReader wig2, String seq_name ) {
		for( BedRecord bedRecord : bedRecords ) {
			if( position >= bedRecord.getStart() && position <= bedRecord.getEnd() ) {
				return true;
			}
			if( bedRecord.getStart() > position ) {
				break;
			}
		}	
		Contig result;
		double thisMean;
		try {
//			
//			result = wig.query(chr, position, position);
//			thisMean = result.mean();
//			if (thisMean > 0){
//				return true;
//			}
			
			result = wig2.query(seq_name, position2, position2);
			thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
		} catch (WigFileException | IOException e) {
			e.printStackTrace();
		}		
		return false;
	}
	
	static boolean ifCallAble( int position, ArrayList<BedRecord> bedRecords, WigFileReader wig, String chr ) {
		for( BedRecord bedRecord : bedRecords ) {
			if( position >= bedRecord.getStart() && position <= bedRecord.getEnd() ) {
				return true;
			}
			if( bedRecord.getStart() > position ) {
				break;
			}
		}	
		Contig result;
		try {
			result = wig.query(chr, position, position);
			double thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
		} catch (WigFileException | IOException e) {
			e.printStackTrace();
		}		
		return false;
	}
}
