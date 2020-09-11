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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;

public class OverlapCnsWithDomestication {
	public static void main(String[] args) {
		try {
			for( int chro=1; chro<=10; ++chro ) {
				String chr = Integer.toString(chro);
				
				HashSet<Integer> allGoodPositions = new HashSet<Integer>();
				HashSet<Integer> cnsGoodPositions = new HashSet<Integer>();
				HashSet<Integer> selectionGoodPositions = new HashSet<Integer>();
				
				TabixReader tr = new TabixReader("/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/B73_callable_status.bed.gz");
				HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46_gene.fa").getChromoSomeHashMap();
				
				
				{
					String s;
					while ((s = tr.readLine()) != null) {
	//					System.out.println(s);
						String[] arrOfStr = s.split("\\s+"); 
						if( arrOfStr[3].compareTo("CALLABLE")==0 && arrOfStr[0].compareTo(chr) ==0 ) {
							for( int posi = Integer.parseInt(arrOfStr[1])+1; posi <=Integer.parseInt(arrOfStr[2]); ++posi  ) {
								//if( maizeGenome.get(chr).charAt(posi-1) != 'd' && maizeGenome.get(chr).charAt(posi-1) != 'D' && maizeGenome.get(chr).charAt(posi-1) != 'N' && maizeGenome.get(chr).charAt(posi-1) != 'n' ) {
									allGoodPositions.add(posi);
								//}
							}
						}
					}
				}
				
				System.out.println("allGoodPositions done ");
				
				SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
				SAMRecordIterator it = reader0.iterator();
				while (it.hasNext()) {
		            final SAMRecord samRecord = it.next();
		            int start = samRecord.getAlignmentStart();
		            int end = samRecord.getAlignmentEnd();
		            if( samRecord.getContig().compareTo(chr)==0 ) {
		            	for( int i=start; i<=end; ++i ) {
			            	if( allGoodPositions.contains(i) ) {
			            		cnsGoodPositions.add(i);
			            	}
			            }
		            }
				}
				reader0.close();
				System.out.println("cnsGoodPositions done ");
				
				
				File file = new File("/media/bs674/1_8t/AndCns/overLapWithDomestication/domesticationLoci/selectionRegions_us_uniq_v4");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = null;
				Pattern p = Pattern.compile("^(\\S*)\t(\\d+)\t(\\d+)\t");
				while ((tempString = reader.readLine()) != null) {
					Matcher m=p.matcher(tempString);
					if(m.find()){
						int start = Integer.parseInt(m.group(2));
	                	int end = Integer.parseInt(m.group(3));
	                	if(  m.group(1).compareTo(chr)==0 ) {
	    	            	for( int i=start+1; i<=end; ++i ) {
	    		            	if( allGoodPositions.contains(i) ) {
	    		            		selectionGoodPositions.add(i);
	    		            	}
	    		            }
	    	            }
					}
				}
				System.out.println("selectionGoodPositions done ");
				reader.close();
				tr.close();
				
				int allLength = 0;
				int cnsLength = 0;
				int selectionLength = 0;
				int cnsLengthAndSelectionLength = 0;
				allLength = allLength + allGoodPositions.size();
				cnsLength = cnsLength + cnsGoodPositions.size();
				selectionLength = selectionLength + selectionGoodPositions.size();
				HashSet<Integer> intersection = new HashSet<Integer>(cnsGoodPositions); // use the copy constructor
				intersection.retainAll(selectionGoodPositions);
				cnsLengthAndSelectionLength = cnsLengthAndSelectionLength + intersection.size();

				System.out.println("allLength:" + allLength + " cnsLength:" + cnsLength + " selectionLength:" + selectionLength + " cnsLengthAndSelectionLength:" + cnsLengthAndSelectionLength);
			}
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
