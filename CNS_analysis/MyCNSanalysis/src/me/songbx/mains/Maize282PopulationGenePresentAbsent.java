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
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.BedRecord;
import me.songbx.model.GeneSimple;
import me.songbx.service.IfIntron;
import htsjdk.tribble.readers.TabixReader;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileReader;

import java.nio.file.Path;
import java.nio.file.Paths;

public class Maize282PopulationGenePresentAbsent {
	public static void main1(String[] args) {
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("./Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("./Zea_mays.AGPv4.34.gff3").getGenes();
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("./gene_sequence.fa");
			for( String geneid : maizeGenes.keySet() ) {
				String seq_name = maizeGenes.get(geneid).getChromeSomeName()+":"+maizeGenes.get(geneid).getStart()+"-"+maizeGenes.get(geneid).getEnd();
				outPut.println(">"+geneid+"\t"+seq_name);
				outPut.println(maizeGenome.get(maizeGenes.get(geneid).getChromeSomeName()).substring(maizeGenes.get(geneid).getStart()-1, maizeGenes.get(geneid).getEnd()));
			}
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter(args[0] + "_genotype");
			//HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			HashMap<String, WigFileReader> wigFileReaders = new HashMap<String, WigFileReader>();
			HashMap<String, WigFileReader> wigFileReaders2 = new HashMap<String, WigFileReader>();
			IfIntron ifIntron = new IfIntron("./Zea_mays.AGPv4.34.gff3");
			HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("./Zea_mays.AGPv4.34.gff3").getGenes();
			outPut.print("geneName\tGeneRange\tlength");
			HashMap<String, String> seqNameToId = new HashMap<String, String>();
	        try {
	        	//reading gene sequence begin
	        	File fastaFile = new File("gene_sequence.fa");
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
            		//String bwFilePosition =   "./mapreadsToGene/" + baseName  + ".bw";
            		String bwFilePosition =   "/public/home/xpsun/521bwabam/geneSeqWig/bigwig/" + baseName  + ".bw";
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
	        for( String geneid: maizeGenes.keySet() ) {
	            int start = maizeGenes.get(geneid).getStart();
	            int end = maizeGenes.get(geneid).getEnd();
	            String chr = maizeGenes.get(geneid).getChromeSomeName();
	            HashSet<Integer> nonMaskingPositions = new HashSet<Integer>();
	        
	            int totalLength = 0;
	            for( int i=start-1; i<end; ++i ) {
	            	if (ifIntron.getElement ( chr, i+1 ) == 2 ) { 
	            		totalLength = totalLength + 1;
        				nonMaskingPositions.add(i);
	            	}
	            }
	            if(totalLength > 30) {
		            outPut.print(geneid+"\t"+chr+":"+start+"-"+end+"\t" + totalLength);
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
		    		            	if( maizeGenome.get(chr).charAt(i) == 'n' && ifCallAble(i+1, bedRecords, wig, chr, i+2-start, wig2, seqNameToId.get(seq_name)) ) {
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
