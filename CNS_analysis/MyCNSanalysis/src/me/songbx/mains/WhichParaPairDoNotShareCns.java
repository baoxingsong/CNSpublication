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

import org.apache.commons.lang3.StringUtils;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.GeneSimple;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class WhichParaPairDoNotShareCns {
	public static void main(String[] args) {
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, GeneSimple> sorghumGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3").getGenes();
		HashMap<String, String> maize2 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
		// get the syntenic paralogous relationship between maize and sorghum begin
		HashMap<String, ArrayList<String>> sorghum_gene_to_maize_genes = new HashMap<String, ArrayList<String>>();
		try {
        	File file = new File("/media/bs674/1_8t/core_Andropogoneae_genes/tabasco2.0/quota-alignment/maize_sorghum/syntenic_genes");
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
			
			while ((tempString = reader.readLine()) != null) {
				if(tempString.startsWith("#")){
					
				}else{
					String[] arrOfStr = tempString.split("\\s+");
					if ( ! sorghum_gene_to_maize_genes.containsKey(arrOfStr[0]) ) {
						sorghum_gene_to_maize_genes.put(arrOfStr[0], new ArrayList<String>());
					}
					sorghum_gene_to_maize_genes.get(arrOfStr[0]).add(arrOfStr[1]);
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		// get the syntenic paralogous relationship between maize and sorghum end
		
		
		PrintWriter outPut = null;
		try {
			int sorghumExtendLength = 2000;
			int maizeExtendLength = 2000;
			int miniscore= 0;
			int maxscore= 50000;
			
			
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/subfunctionlization/sorghum/subfunctionlization_V2_not_share_CNS"+miniscore+"_"+maxscore+"_"+sorghumExtendLength+"_"+maizeExtendLength);
	        Pattern p = Pattern.compile("^(\\d+)H");
	        Pattern pattern = Pattern.compile("(\\d+)(\\w)");
			for ( String sorghumGene : sorghum_gene_to_maize_genes.keySet() ) {
				if( sorghum_gene_to_maize_genes.get(sorghumGene).size() == 2 ) {
					{
						int i=0;
						if ( ! sorghum_gene_to_maize_genes.containsKey(sorghumGene) ) {
							System.out.println("line 63:" + sorghumGene);
							continue;
						}
						if ( ! maizeGenes.containsKey(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)) ) {
							System.out.println("line 67:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
						if ( ! sorghumGenes.containsKey(sorghumGene) ) {
							System.out.println("line 71:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
						
						i=1;
						if ( ! maizeGenes.containsKey(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)) ) {
							System.out.println("line 82:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
						if ( ! sorghumGenes.containsKey(sorghumGene) ) {
							System.out.println("line 85:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
					}
					
					{ // all the cis region
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/2t/testPan_cns/testAlgorithm/sorghum.bam"));
						SAMRecordIterator samit = reader0.queryOverlapping(sorghumIntervalChr, sorghumGenes.get(sorghumGene).getStart(), sorghumGenes.get(sorghumGene).getEnd());
						
						boolean find0 = false;
						boolean find1 = false;
						
						while (samit.hasNext()) {
							final SAMRecord samRecord = samit.next();
							if( samRecord.getReadName().compareTo(sorghum_gene_to_maize_genes.get(sorghumGene).get(0)) == 0 ) {
								find0 = true;
							}
							if( samRecord.getReadName().compareTo(sorghum_gene_to_maize_genes.get(sorghumGene).get(1)) == 0 ) {
								find1 = true;
							}
						}
						reader0.close();
						if ( find0 && find1 ) {
							//System.out.println("line 110");
							int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getStart() - sorghumExtendLength;
							int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getEnd() + sorghumExtendLength;
							
							HashSet<Integer> allPositions = new HashSet<Integer>();
							HashSet<Integer> positions1 = new HashSet<Integer>();
							HashSet<Integer> positions2 = new HashSet<Integer>();
							
							for (int i=0; i<2; ++i) {
								int wantedStrand = 1;
								if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
									wantedStrand = 0;
								}
								
								String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
								int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeExtendLength;
								int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd() + maizeExtendLength;
								
								SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
								
								SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
								while (it.hasNext()) {
						            final SAMRecord samRecord = it.next();
						            Matcher m=p.matcher(samRecord.getCigarString());
						            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore ){
						            	int strand = 1;
							            if ( samRecord.getReadNegativeStrandFlag() ) {
							            	strand = 0;
							            }
							            if ( strand == wantedStrand ) {
							            	int queryStart = Integer.parseInt(m.group(1));
							            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
							            	if( samRecord.getReadName().compareTo(sorghumIntervalChr)==0 ) {
							            		if( (sorghumIntervalStart<=queryStart && queryStart<= sorghumIntervalEnd) ||
							            				(sorghumIntervalStart<=queryEnd && queryEnd<= sorghumIntervalEnd) ||
							            				(queryStart<=sorghumIntervalEnd && sorghumIntervalEnd<= queryEnd) ||
							            				(queryStart<=sorghumIntervalStart && sorghumIntervalStart<= queryEnd)
							            				) {
							            			String seqMaize = maize2.get(maizeChr).substring(samRecord.getAlignmentStart()-1, samRecord.getAlignmentEnd());
							            			String seqsorghum = samRecord.getReadString();
							            			String alnMaize = "";
							            			String alnSorghum = "";
							            			int maizePosition = 0;
							            			int sorghumPostion = 0;
							            			Matcher matcher = pattern.matcher(samRecord.getCigarString());
									            	while( matcher.find() ) {
									            		if( matcher.group(2).compareTo("H") == 0 ) {
									            			
									            		}else if( matcher.group(2).compareTo("M") == 0 ) {
									            			int length = Integer.parseInt(matcher.group(1));
									            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
									            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
									            			maizePosition += length;
									            			sorghumPostion += length;
									            		}else if( matcher.group(2).compareTo("I") == 0 ) {
									            			int length = Integer.parseInt(matcher.group(1));
									            			alnMaize += StringUtils.repeat("-", length);
									            			alnSorghum += seqsorghum.substring(sorghumPostion, sorghumPostion+length);
									            			sorghumPostion += length;
									            		}else if( matcher.group(2).compareTo("D") == 0 ) {
									            			int length = Integer.parseInt(matcher.group(1));
									            			alnMaize += seqMaize.substring(maizePosition, maizePosition+length);
									            			alnSorghum += StringUtils.repeat("-", length);
									            			maizePosition += length;
									            		}else {
									            			System.out.println("some thing unknow happened");
									            		}
									            	}
									            	if(samRecord.getFlags()==0) {
									            		int position = queryStart;
									            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
															if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip) && position>=sorghumIntervalStart && position<=sorghumIntervalEnd )  {
																if( 0 == i ) {
																	positions1.add(position);
																}else {
																	positions2.add(position);
																}
																allPositions.add(position);
															}
															if( alnSorghum.charAt(ip)  != '-') {
																++position;
															}
														}
									            	}else{
									            		int position = queryEnd;
									            		for( int ip = 0; ip<alnMaize.length(); ++ip ) {
															if( alnMaize.charAt(ip) != 'N' && alnMaize.charAt(ip)==alnSorghum.charAt(ip) && position>=sorghumIntervalStart && position<=sorghumIntervalEnd )  {
																if( 0 == i ) {
																	positions1.add(position);
																}else {
																	positions2.add(position);
																}
																allPositions.add(position);
															}
															if( alnSorghum.charAt(ip)  != '-') {
																--position;
															}
														}
									            	}
							            		}
							            	}
							            }
									}
						        }
								try {
									reader.close();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
							if( (positions1.size() == 0 || positions2.size()==0) && ( positions1.size() != 0 || positions2.size() != 0 ) ) {
								outPut.println("all\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + allPositions.size() + "\t" + positions1.size() + "\t" + positions2.size());
							}else if ( (positions1.size() == 0 || positions2.size()==0)  ) {
								System.out.println( "" + positions1.size() + " " + positions2.size());
							}
						}
					}
				}
			}
			outPut.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
	}
}
