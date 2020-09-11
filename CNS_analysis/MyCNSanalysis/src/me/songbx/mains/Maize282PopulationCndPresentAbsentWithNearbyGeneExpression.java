package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

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


public class Maize282PopulationCndPresentAbsentWithNearbyGeneExpression {
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSBasedGwas/genotype");
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			HashMap<String, WigFileReader> wigFileReaders = new HashMap<String, WigFileReader>();
			
			outPut.print("GenetypeID\tlength");
	        try {
	        	File file = new File("/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/callable_list");
	        	BufferedReader reader = new BufferedReader(new FileReader(file));
	            String tempString = null;
	            while ((tempString = reader.readLine()) != null) {
	            	tempString = tempString.trim();
	        		String baseName = tempString;
	        		baseName = baseName.replaceAll(".*\\/", "");
	        		baseName = baseName.replaceAll("_callable_status.bed.gz", "");
	        		outPut.print("\t" + baseName);
	        		
	        		String bwFilePosition =   "/media/bs674/panAndAssemblyfi/maize282Genotyping/coverageBwFiles/" + baseName  + ".bw";
            		Path bwFile = Paths.get(bwFilePosition);
            		WigFileReader wig = WigFileReader.autodetect(bwFile);
            		wigFileReaders.put(baseName, wig);
            		
					TabixReader tr = new TabixReader(tempString);
					tabixReaders.put(baseName, tr);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        outPut.println();
	        
	        String bwFilePosition = "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw";
	        Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig0 = WigFileReader.autodetect(bwFile);
			
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/Setaria_italica/result/5.bam"));
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
	            outPut.print(chr+":"+start+"-"+end+"\t" + totalLength);
	            //read bed files begin
	    		File file = new File("/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/callable_list");
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
			try {
				reader0.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			for( String baseName : tabixReaders.keySet()  ) {
				tabixReaders.get(baseName).close();
				wigFileReaders.get(baseName).close();
			}
			wig0.close();
		}catch (IOException e) {
            e.printStackTrace();
        } catch (WigFileException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
