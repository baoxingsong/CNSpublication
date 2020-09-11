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
import me.songbx.model.Strand;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CnsSubfunctionlizationV2 {
	public static void main(String[] args) {
		new CnsSubfunctionlizationV2(1500);
		new CnsSubfunctionlizationV2(1600);
		new CnsSubfunctionlizationV2(1700);
		new CnsSubfunctionlizationV2(1800);
		new CnsSubfunctionlizationV2(1900);
		new CnsSubfunctionlizationV2(2000);
		new CnsSubfunctionlizationV2(2100);
		new CnsSubfunctionlizationV2(2200);
		new CnsSubfunctionlizationV2(2300);
		new CnsSubfunctionlizationV2(2400);
		new CnsSubfunctionlizationV2(2500);
	}
	public CnsSubfunctionlizationV2( int maizeIntervalEnd){
		int maizeIntervalStart = 1;
		
		int sorghumExtendLength = maizeIntervalEnd;
		int miniscore = 0;
		int maxscore = 500000;
		
		// read gff files begin
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, GeneSimple> sorghumGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3").getGenes();
		// read gff files end
		
		//read genome files begin
//		HashMap<String, String> sorghum1 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_sorghum_k20_26.fa").getChromoSomeHashMap();
		HashMap<String, String> sorghum2 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_sorghum_k20_26_gff_cds.fa").getChromoSomeHashMap();
		//HashMap<String, String> maize = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
		// read genome files end
		
		// generated a data set of successfully lifted maize genes begin
		HashMap<String, HashSet<String>> liftedMaizeToSorghum = new HashMap<String, HashSet<String>>(); //key is the maize gene and value is a hashset of sorghum genes
		try {
			SamReader readerSam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/sorghum.sam"));
			SAMRecordIterator itSam = readerSam.iterator();
			while (itSam.hasNext()) {
	            final SAMRecord samRecord = itSam.next();
	            if( samRecord.getFlags()!= 4) {
	            	if( ! liftedMaizeToSorghum.containsKey(samRecord.getReadName()) ) {
		            	liftedMaizeToSorghum.put(samRecord.getReadName(), new HashSet<String>());
		            }
		            String sorghumChr = samRecord.getContig();
		            int start = samRecord.getStart();
		            int end = samRecord.getEnd();
//		            System.out.println(samRecord.getFlags() + " \t" + sorghumChr);
		            for( String geneId : sorghumGenes.keySet()  ) {
		            	if ( sorghumGenes.get(geneId).getChromeSomeName().compareTo(sorghumChr) == 0 ) {
		            		if( (sorghumGenes.get(geneId).getStart() <= start && start <= sorghumGenes.get(geneId).getEnd()) || (
		            				sorghumGenes.get(geneId).getStart() <= end && end <= sorghumGenes.get(geneId).getEnd() )) {
		            			if( (sorghumGenes.get(geneId).getStrand() == Strand.NEGTIVE && samRecord.getReadNegativeStrandFlag() ) ||
		            					(sorghumGenes.get(geneId).getStrand() != Strand.NEGTIVE && (!samRecord.getReadNegativeStrandFlag() ) )	) {
		            				liftedMaizeToSorghum.get(samRecord.getReadName()).add(geneId);
//		            				System.out.println(samRecord.getReadName() + "\t" + geneId);
		            			}
		            		}
		            	}
		            }
	            }
			}
			readerSam.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		// generated a data set of successfully lifted maize genes begin
		
		
		// get the syntenic paralogous relationship between maize and sorghum begin
		HashMap<String, ArrayList<String>> sorghum_gene_to_maize_genes = new HashMap<String, ArrayList<String>>();
		try {
        	File file = new File("/media/bs674/pan_and_non_asse/core_Andropogoneae_genes/tabasco2.0/quota-alignment/maize_sorghum/syntenic_genes");
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
			
			while ((tempString = reader.readLine()) != null) {
				if(tempString.startsWith("#")){
					
				}else{
					String[] arrOfStr = tempString.split("\\s+");
					if ( ! sorghum_gene_to_maize_genes.containsKey(arrOfStr[0]) ) {
						sorghum_gene_to_maize_genes.put(arrOfStr[0], new ArrayList<String>());
					}
					if( liftedMaizeToSorghum.containsKey(arrOfStr[1]) && liftedMaizeToSorghum.get(arrOfStr[1]).contains(arrOfStr[0]) ) {
						sorghum_gene_to_maize_genes.get(arrOfStr[0]).add(arrOfStr[1]);
					}
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		
		
		
		PrintWriter outPut = null;
		try {
			/*
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			
			try {
				File file = new File( "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.wig" );
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = null;
				String chr = "";
				Pattern p0 =Pattern.compile("chrom=(\\S+)");
				while ((tempString = reader.readLine()) != null) {
					if( tempString.charAt(0) == 't' ) {
						
					}else if( tempString.charAt(0) == 'v' ) {
						Matcher m = p0.matcher(tempString);
						if( m.find() ) {
							chr = m.group(1);
						}
					}else {
						String[] arrOfStr = tempString.split("\\s+");
						if( ! cnsRecords.containsKey(chr) ) {
							cnsRecords.put(chr, new HashSet<Integer>());
						}
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							cnsRecords.get(chr).add(Integer.parseInt(arrOfStr[0]));
						}
						
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			*/
			
			
			
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/subfunctionlization/pancns/subfunctionlization_score"+miniscore+"_"+maxscore+"_"+sorghumExtendLength+"_"+maizeIntervalStart+"_"+maizeIntervalEnd);
	        Pattern p = Pattern.compile("^(\\d+)H");
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

					
					
					// this also count number of Ms with cigar string , similar with samtools depth command
					int all_position=0;
					int all_position1=0;
					int all_position2=0;
					{ // up stream, upstream is defined from start codon
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getStart() - sorghumExtendLength;
						int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getStart()-1;
						
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment2>();
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment2>();
						for (int i=0; i<2; ++i) {
//							System.out.println("iiiii:" + i);
							int wantedStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
								wantedStrand = 0; // negative strand
							}
							
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeIntervalEnd;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-maizeIntervalStart;
							if ( 0 == wantedStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalStart;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalEnd;
							}
							SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
							
							SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
							while (it.hasNext()) {
					            final SAMRecord samRecord = it.next();
					            
//					            System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
			            		
					            Matcher m=p.matcher(samRecord.getCigarString());
					            if(m.find() && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore){
					            	int strand = 1;
						            if ( samRecord.getReadNegativeStrandFlag() ) {
						            	strand = 0;
						            }
						            if ( strand == wantedStrand ) {
//						            	System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
					            		
						            	int queryStart = Integer.parseInt(m.group(1));
						            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
//						            	System.out.println("" + samRecord.getReadString().length() + "\t" + queryStart + " " + queryEnd);
					            		
						            	if( samRecord.getReadName().compareTo(sorghumIntervalChr)==0 && sorghumIntervalChr.charAt(0)!='s' &&  maizeChr.charAt(0)!='B' ) {
						            		if( (sorghumIntervalStart<=queryStart && queryStart<= sorghumIntervalEnd) ||
						            				(sorghumIntervalStart<=queryEnd && queryEnd<= sorghumIntervalEnd) ||
						            				(queryStart<=sorghumIntervalEnd && sorghumIntervalEnd<= queryEnd) ||
						            				(queryStart<=sorghumIntervalStart && sorghumIntervalStart<= queryEnd) ) {
						            			//System.out.println( maizeChr + "\t" + samRecord.getStart() + "\t" + sorghumIntervalChr + "\t" + queryStart + "\t" + samRecord.getCigarString() );
//						            			ArrayList<String> alignment = cirgarToAlignment( maize.get(maizeChr), samRecord.getStart(), sorghum2.get(sorghumIntervalChr), queryStart, samRecord.getCigarString()  );
//						            			System.out.println(alignment.get(0));
//						            			System.out.println(alignment.get(1) + "\n");
						            			AlignedSorghumFragment2 alignedSorghumFragment = new AlignedSorghumFragment2(samRecord.getReadName(), strand);
						            			if ( samRecord.getReadNegativeStrandFlag() ) {
						            				for ( int position = queryEnd; position >= queryStart; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
							            								alignedSorghumFragment.getPositions().add(position);
											            			}
							            							--position;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						position-=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}else {
						            				for ( int position = queryStart; position<=queryEnd; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
							            								alignedSorghumFragment.getPositions().add(position);
											            			}
							            							++position;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						position+=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}
						            			
								            	if( 0 == i ) {
								            		alignedSorghumFragments1.add(alignedSorghumFragment);
								            	}else {
								            		alignedSorghumFragments2.add(alignedSorghumFragment);
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
						HashSet<Integer> allPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions2.add(position);
								}
							}
						}
						if ( Strand.POSITIVE == sorghumGenes.get(sorghumGene).getStrand()  ) {
							outPut.print("upstream");
						}else {
							outPut.print("downstream");
						}
						all_position += allPositions.size();
						all_position1 += positions1.size();
						all_position2 += positions2.size();
						outPut.println("\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + allPositions.size() + "\t" + positions1.size() + "\t" + positions2.size());
					}
					
					if( sorghumGenes.get(sorghumGene).getNumberOfCds() > 1 ) { // intron
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getStart();
						int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getEnd();
						
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment2>();
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment2>();
						for (int i=0; i<2; ++i) {
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getNumberOfCds() > 1 ) {
								int wantedStrand = 1;
								if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
									wantedStrand = 0;
								}
								
								String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
								int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart();
								int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd();
								SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
								
								SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
								while (it.hasNext()) {
						            final SAMRecord samRecord = it.next();
						            
						            Matcher m=p.matcher(samRecord.getCigarString());
						            if(m.find()  && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore){
						            	int strand = 1;
							            if ( samRecord.getReadNegativeStrandFlag() ) {
							            	strand = 0;
							            }
							            if ( strand == wantedStrand ) {
							            	//System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
								            
							            	int queryStart = Integer.parseInt(m.group(1));
							            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
							            	if( samRecord.getReadName().compareTo(sorghumIntervalChr)==0 ) {
							            		//System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
							            		//System.out.println("" + samRecord.getReadString().length() + "\t" + queryStart + " " + queryEnd);
							            		if( (sorghumIntervalStart<=queryStart && queryStart<= sorghumIntervalEnd) ||
							            				(sorghumIntervalStart<=queryEnd && queryEnd<= sorghumIntervalEnd) ||
							            				(queryStart<=sorghumIntervalEnd && sorghumIntervalEnd<= queryEnd) ||
							            				(queryStart<=sorghumIntervalStart && sorghumIntervalStart<= queryEnd)
							            				
							            				) {
							            			AlignedSorghumFragment2 alignedSorghumFragment = new AlignedSorghumFragment2(samRecord.getReadName(), strand);
							            			if ( samRecord.getReadNegativeStrandFlag() ) {
							            				for ( int position = queryEnd; position >= queryStart; ) {
								            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
								            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
								            				while (matcher.find()) {
								            					int length = Integer.parseInt( matcher.group(1));
								            					char cigerLetter =  matcher.group(2).charAt(0);
								            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
								            						for( int il=0; il<length; ++il ) {
								            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
								            								alignedSorghumFragment.getPositions().add(position);
												            			}
								            							--position;
								            						}
								            					} else if ( cigerLetter=='I' ) {
								            						position-=length;
								            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
								            						
								            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
								            					
								            					}else {
								            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
								            					}
								            				}
								            			}
							            			}else {
							            				for ( int position = queryStart; position<=queryEnd; ) {
								            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
								            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
								            				while (matcher.find()) {
								            					int length = Integer.parseInt( matcher.group(1));
								            					char cigerLetter =  matcher.group(2).charAt(0);
								            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
								            						for( int il=0; il<length; ++il ) {
								            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
								            								alignedSorghumFragment.getPositions().add(position);
												            			}
								            							++position;
								            						}
								            					} else if ( cigerLetter=='I' ) {
								            						position+=length;
								            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
								            						
								            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
								            					
								            					}else {
								            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
								            					}
								            				}
								            			}
							            			}
									            	if( 0 == i ) {
									            		alignedSorghumFragments1.add(alignedSorghumFragment);
									            	}else {
									            		alignedSorghumFragments2.add(alignedSorghumFragment);
									            	}
							            		}
							            	}
							            }
									}
						        }
								//System.out.println();
								//System.out.println();
								try {
									reader.close();
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}
						}
						HashSet<Integer> allPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions2.add(position);
								}
							}
						}
						outPut.println("intron\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + allPositions.size() + "\t" + positions1.size() + "\t" + positions2.size());
					}
					{ // down stream
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getEnd()+1;
						int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getEnd()+sorghumExtendLength;
						
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment2>();
						ArrayList<AlignedSorghumFragment2> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment2>();
						for (int i=0; i<2; ++i) {
							int wantedStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
								wantedStrand = 0;
							}
							
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalStart;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalEnd;
							
							if ( 0 == wantedStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeIntervalEnd;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-maizeIntervalStart;
							}
							
							SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
							SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
							while (it.hasNext()) {
					            final SAMRecord samRecord = it.next();
					            Matcher m=p.matcher(samRecord.getCigarString());
					            if(m.find()  && samRecord.getMappingQuality() >= miniscore && samRecord.getMappingQuality()<=maxscore){
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
						            			AlignedSorghumFragment2 alignedSorghumFragment = new AlignedSorghumFragment2(samRecord.getReadName(), strand);
						            			if ( samRecord.getReadNegativeStrandFlag() ) {
						            				for ( int position = queryEnd; position >= queryStart; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
							            								alignedSorghumFragment.getPositions().add(position);
											            			}
							            							--position;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						position-=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}else {
						            				for ( int position = queryStart; position<=queryEnd; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( position >= sorghumIntervalStart && position<=sorghumIntervalEnd ) {
							            								alignedSorghumFragment.getPositions().add(position);
											            			}
							            							++position;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						position+=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}
								            	if( 0 == i ) {
								            		alignedSorghumFragments1.add(alignedSorghumFragment);
								            	}else {
								            		alignedSorghumFragments2.add(alignedSorghumFragment);
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
						HashSet<Integer> allPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment2 alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.getPositions() ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n' && sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'N')  {
									allPositions.add(position);
									positions2.add(position);
								}
							}
						}
						if ( Strand.POSITIVE == sorghumGenes.get(sorghumGene).getStrand()  ) {
							outPut.print("downstream");
						}else {
							outPut.print("upstream");
						}
						all_position += allPositions.size();
						all_position1 += positions1.size();
						all_position2 += positions2.size();
						outPut.println("\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + allPositions.size() + "\t" + positions1.size() + "\t" + positions2.size());
					}
					outPut.println("all\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + all_position + "\t" + all_position1 + "\t" + all_position2);
				}
			}
			outPut.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
	}
	public static ArrayList<String> cirgarToAlignment( String referenceSequence, int referenceBpPos, String querySequence, int queryBpPos, String cigar  ){
		String queryAlignment="";
		String referenceAlignment="";
		Pattern p = Pattern.compile("(\\d+)([MIDNSHP=XB])");
		Matcher matcher = p.matcher(cigar);
		while (matcher.find()) {
			int length = Integer.parseInt( matcher.group(1));
			char cigerLetter =  matcher.group(2).charAt(0);
			if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
				for( int il=0; il<length; ++il ) {
					queryAlignment += querySequence.substring(queryBpPos-1, queryBpPos);
					referenceAlignment += referenceSequence.substring(referenceBpPos-1, referenceBpPos);
					++queryBpPos;
					++referenceBpPos;
				}
			} else if ( cigerLetter=='I' ) {
				referenceAlignment += StringUtils.repeat("-", length);
				for( int il=0; il<length; ++il ) {
					referenceAlignment += StringUtils.repeat("-", length);
					queryAlignment += querySequence.substring(queryBpPos-1, queryBpPos);
					++queryBpPos;
				}
			} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
				referenceAlignment += referenceSequence.substring(referenceBpPos-1, referenceBpPos+length-1);
				queryAlignment += StringUtils.repeat("-", length);
				referenceBpPos+=length;
			} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
			}else {
				System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
			}
		}
		ArrayList<String> alignment = new ArrayList<String>();
		alignment.add(referenceAlignment);
		alignment.add(queryAlignment);
		return alignment;
	}
}

class AlignedSorghumFragment2{
	private String chr;
	private HashSet<Integer> positions;
	private int strand; //1 positive, 0 negative
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public HashSet<Integer> getPositions() {
		return positions;
	}
	public void setPositions(HashSet<Integer> positions) {
		this.positions = positions;
	}
	public int getStrand() {
		return strand;
	}
	public void setStrand(int strand) {
		this.strand = strand;
	}
	public AlignedSorghumFragment2(String chr, int strand) {
		super();
		this.chr = chr;
		this.positions = new HashSet<Integer>();
		this.strand = strand;
	}
}
