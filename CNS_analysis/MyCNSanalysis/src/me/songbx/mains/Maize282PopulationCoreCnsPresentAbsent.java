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
import edu.unc.genomics.io.WigFileReader;

import java.nio.file.Path;
import java.nio.file.Paths;

public class Maize282PopulationCoreCnsPresentAbsent {
	
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter(args[0] + "_genotype");
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			HashMap<String, WigFileReader> wigFileReaders = new HashMap<String, WigFileReader>();
			HashMap<String, WigFileReader> wigFileReaders2 = new HashMap<String, WigFileReader>();

			outPut.print("cns_id\tGeneRange\tlength");
			HashMap<String, String> seqNameToId = new HashMap<String, String>();
	        try {
	        	//reading gene sequence begin
	        	File fastaFile = new File("cns_seq.fa");
	        	BufferedReader fastaReader = new BufferedReader(new FileReader(fastaFile));
	        	String tempString = null;
	        	Pattern pattern = Pattern.compile("^>(\\S+)\\s+(\\S+)$");
	        	while ((tempString = fastaReader.readLine()) != null) {
	        		Matcher match = pattern.matcher(tempString);
	            	if( match.find() ) {
	            		seqNameToId.put(match.group(2), match.group(1));
	            	}
	        	}
	        	fastaReader.close();
	        	//reading gene sequence done
	        	
	        	//reading reads coverage begin
	        	File file = new File(args[0]);
	        	BufferedReader reader = new BufferedReader(new FileReader(file));
	            tempString = null;
	            while ((tempString = reader.readLine()) != null) {
	            	tempString = tempString.trim();
	        		String baseName = tempString;
	        		baseName = baseName.replaceAll(".*\\/", "");
	        		baseName = baseName.replaceAll(".bam.bed.gz", "");
	        		outPut.print("\t" + baseName);
	        		
	        		//String bwFilePosition0 =   "/media/bs674/panAndAssemblyfi/maize282Genotyping/coverageBwFiles/" + baseName  + ".bw";
	        		//String bwFilePosition0 =   "./coverageBwFiles/" + baseName  + ".bw";
	        		String bwFilePosition0 = "/public/home/xpsun/521bwabam/bamtoWig/bigwig/" + baseName + ".bam.wig.bw";
            		Path bwFile0 = Paths.get(bwFilePosition0);
            		WigFileReader wig0 = WigFileReader.autodetect(bwFile0);
            		wigFileReaders.put(baseName, wig0);
            		
            		//String bwFilePosition =   "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/mapreadsToCns/" + baseName  + ".bw";
            		//String bwFilePosition =   "./mapreadsToCns/" + baseName  + ".bw";
            		String bwFilePosition = "/public/home/xpsun/521bwabam/cnswig/wig2big/" + baseName + ".bam.wig.bw";
            		Path bwFile = Paths.get(bwFilePosition);
            		WigFileReader wig = WigFileReader.autodetect(bwFile);
            		wigFileReaders2.put(baseName, wig);
            		
					TabixReader tr = new TabixReader(tempString);
					tabixReaders.put(baseName, tr);
	            }
	            reader.close();
	          //reading reads coverage done
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        outPut.println();
	        
	        //SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
	        for( String seqName: seqNameToId.keySet() ) {
	        	Pattern p = Pattern.compile("^(\\w+):(\\d+)-(\\d+)$");
	        	Matcher m = p.matcher(seqName);
	        	if ( m.find() ) {
		            int start = Integer.parseInt(m.group(2))+1;
		            int end = Integer.parseInt(m.group(3));
		            String chr = m.group(1);
		            
		            HashSet<Integer> nonMaskingPositions = new HashSet<Integer>();
		            int totalLength = 0;
		            for( int i=start-1; i<end; ++i ) {
		            	if ( maizeGenome.get(chr).charAt(i) != 'n' && maizeGenome.get(chr).charAt(i) != 'N' ) {
		            		nonMaskingPositions.add(i);
		            		totalLength++;
		            	}
		            }
		            if ( totalLength > 6) {
			            outPut.print(chr+":"+start+"-"+end+"\t" + totalLength);
		            
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
			                		baseName = baseName.replaceAll(".bam.bed.gz", "");
			                		WigFileReader wig = wigFileReaders.get(baseName);
			                		WigFileReader wig2 = wigFileReaders2.get(baseName);
			                		TabixReader tr = tabixReaders.get(baseName);
		
			                		String s;
			                		int st = start-1;
			    					TabixReader.Iterator iter = tr.query(chr+":"+st+"-"+end); // get the iterator
			    					while ((s = iter.next()) != null) {
			    						String[] arrOfStr = s.split("\\s+"); 
			    						if( arrOfStr[3].compareTo("CALLABLE")==0 ) {
			    							BedRecord bedRecord = new BedRecord(Integer.parseInt(arrOfStr[1])+1, Integer.parseInt(arrOfStr[2])); //bed has special coordinate
			    							bedRecords.add(bedRecord);
			    						}
			    					}
			    					int totalGoodLength = 0;
			    		            for( int i : nonMaskingPositions ) {
			    		            	if( ifCallAble(i+1, bedRecords, wig, chr, i+2-start, wig2, seqNameToId.get(seqName)) ) {
			    		            		totalGoodLength++;
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
        }// catch (WigFileException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
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
			result = wig.query(chr, position, position);
			thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
			
			result = wig2.query(seq_name, position2, position2);
			thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}		
		return false;
	}
}
