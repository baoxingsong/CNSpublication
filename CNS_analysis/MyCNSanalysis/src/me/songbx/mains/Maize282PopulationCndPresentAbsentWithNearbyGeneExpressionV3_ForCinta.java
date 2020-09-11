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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;

import java.nio.file.Path;
import java.nio.file.Paths;


public class Maize282PopulationCndPresentAbsentWithNearbyGeneExpressionV3_ForCinta {
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter(args[0] + "_genotype");
			//HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			HashMap<String, WigFileReader> wigFileReaders = new HashMap<String, WigFileReader>();
			HashMap<String, WigFileReader> wigFileReaders2 = new HashMap<String, WigFileReader>();
			
			outPut.print("GenetypeID\tlength");
			HashMap<String, String> seqNameToId = new HashMap<String, String>();
	        try {
	        	//File fastaFile = new File("/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/cns_seq.fa");
	        	File fastaFile = new File("cns_seq.fa");
	        	BufferedReader fastaReader = new BufferedReader(new FileReader(fastaFile));
	        	String tempString = null;
	        	Pattern pattern = Pattern.compile("^>(\\S+)\\s+(\\S+)");
	        	while ((tempString = fastaReader.readLine()) != null) {
	        		Matcher match = pattern.matcher(tempString);
	            	if( match.find() ) {
	            		seqNameToId.put(match.group(2), match.group(1));
	            		//System.out.println(match.group(2)+"\t"+ match.group(1));
	            	}
	        	}
	        	fastaReader.close();
	        	
	        	File file = new File(args[0]);
	        	BufferedReader reader = new BufferedReader(new FileReader(file));
	            
	            while ((tempString = reader.readLine()) != null) {
	            	tempString = tempString.trim();
	            	String[] currencies = tempString.split("\\s+");
	        		String baseName = currencies[0];
	        		baseName = baseName.replaceAll(".*\\/", "");
	        		baseName = baseName.replaceAll("_status.bed.gz", "");
	        		outPut.print("\t" + baseName);
	        		
	        		//String bwFilePosition0 =   "/media/bs674/panAndAssemblyfi/maize282Genotyping/coverageBwFiles/" + baseName  + ".bw";
	        		String bwFilePosition0 =   "./" + baseName  + ".bw";
            		Path bwFile0 = Paths.get(bwFilePosition0);
            		WigFileReader wig0 = WigFileReader.autodetect(bwFile0);
            		wigFileReaders.put(baseName, wig0);
            		
            		//String bwFilePosition =   "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/mapreadsToCns/" + baseName  + ".bw";
            		String bwFilePosition =   currencies[1];
            		Path bwFile = Paths.get(bwFilePosition);
            		WigFileReader wig = WigFileReader.autodetect(bwFile);
            		wigFileReaders2.put(baseName, wig);
            		
					TabixReader tr = new TabixReader(currencies[0]);
					tabixReaders.put(baseName, tr);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        outPut.println();
	        
	        //String bwFilePosition = "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw";
	        String bwFilePosition = "./5.bw";
	        Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig0 = WigFileReader.autodetect(bwFile);
	        
	        //SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
	        SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("./5.bam"));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            String chr = samRecord.getContig();
	            HashSet<Integer> nonMaskingPositions = new HashSet<Integer>();
	        
	            int totalLength = 0;
	            for( int i=start-1; i<end; ++i ) {
	            	if( maizeGenome.get(chr).charAt(i) != 'n' ) {
	            		Contig result = wig0.query(chr, i+1, i+1);
	        			double thisMean = result.mean();
	        			if( thisMean>0 ) {
	        				totalLength = totalLength + 1;
	        				nonMaskingPositions.add(i);
	        			}
	            	}
	            }
	            if(totalLength > 30) {
		            outPut.print(chr+":"+start+"-"+end+"\t" + totalLength);
		            String seq_name = chr+":"+start+"-"+end;
		            //System.out.println(seq_name);
		            
		            //read bed files begin
		    		File file = new File(args[0]);
		            try {
		            	BufferedReader reader = new BufferedReader(new FileReader(file));
		                String tempString = null;
		                while ((tempString = reader.readLine()) != null) {
		                	tempString = tempString.trim();
		                	try {
		                		ArrayList<BedRecord> bedRecords = new ArrayList<BedRecord>();
		                		tempString = tempString.trim();
		    	            	String[] currencies = tempString.split("\\s+");
		    	        		String baseName = currencies[0];
		                		baseName = baseName.replaceAll(".*\\/", "");
		                		baseName = baseName.replaceAll("_status.bed.gz", "");
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
		    		            	if( maizeGenome.get(chr).charAt(i) != 'n' && ifCallAble(i+1, bedRecords, wig, chr, i+2-start, wig2, seqNameToId.get(seq_name)) ) {
		    		            	//if( maizeGenome.get(chr).charAt(i) != 'n' && ifCallAble(i+1, bedRecords, wig, chr) ) {
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
            outPut.close();
            wig0.close();
            reader0.close();
			for( String baseName : tabixReaders.keySet()  ) {
				tabixReaders.get(baseName).close();
				wigFileReaders.get(baseName).close();
				wigFileReaders2.get(baseName).close();
			}
		}catch (IOException e) {
            e.printStackTrace();
        } catch (WigFileException e) {
			// TODO Auto-generated catch block
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
