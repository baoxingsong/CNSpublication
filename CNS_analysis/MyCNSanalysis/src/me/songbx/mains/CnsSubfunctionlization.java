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
import me.songbx.model.GeneSimple;
import me.songbx.model.Strand;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CnsSubfunctionlization {
	public static void main(String[] args) {
		// read gff files begin
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, GeneSimple> sorghumGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3").getGenes();
		// read gff files end
		
		//read genome files begin
		HashMap<String, String> sorghum1 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_sorghum_k20_26.fa").getChromoSomeHashMap();
		HashMap<String, String> sorghum2 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_sorghum_k20_26_gff_cds.fa").getChromoSomeHashMap();
		// read genome files end
		
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
					sorghum_gene_to_maize_genes.get(arrOfStr[0]).add(arrOfStr[1]);
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		// get the syntenic paralogous relationship between maize and sorghum end
		
		// TODO for test
//		sorghum_gene_to_maize_genes.clear();
//		sorghum_gene_to_maize_genes.put("SORBI_3001G121600", new ArrayList<String>());
//		sorghum_gene_to_maize_genes.get("SORBI_3001G121600").add("Zm00001d033673");
//		sorghum_gene_to_maize_genes.get("SORBI_3001G121600").add("Zm00001d013467");
		
		PrintWriter outPut = null;
		try {
			
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
			
			
			int sorghumExtendLength = 20000;
			int maizeExtendLength = 20000;
			int miniscore= 0;
			int maxscore= 500000;
			
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/subfunctionlization/pancns/subfunctionlization_score"+miniscore+"_"+maxscore+"_"+sorghumExtendLength+"_"+maizeExtendLength);
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

					
					
					
					
					// this also count number of bps with reads covered , similar with samtools mpileup command
					int all_position=0;
					int all_position1=0;
					int all_position2=0;
					{ // up stream, upstream is defined from start codon
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getStart() - sorghumExtendLength;
						int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getStart()-1;
						
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment>();
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment>();
						for (int i=0; i<2; ++i) {
//							System.out.println("iiiii:" + i);
							int wantedStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
								wantedStrand = 0;
							}
							
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeExtendLength;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-1;
							if ( 0 == wantedStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+1;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeExtendLength;
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
					            		
						            	if( samRecord.getReadName().compareTo(sorghumIntervalChr)==0 ) {
						            		if( (sorghumIntervalStart<=queryStart && queryStart<= sorghumIntervalEnd) ||
						            				(sorghumIntervalStart<=queryEnd && queryEnd<= sorghumIntervalEnd) ||
						            				(queryStart<=sorghumIntervalEnd && sorghumIntervalEnd<= queryEnd) ||
						            				(queryStart<=sorghumIntervalStart && sorghumIntervalStart<= queryEnd) ) {
						            			if( queryStart < sorghumIntervalStart ) {
						            				queryStart = sorghumIntervalStart;
						            			}
						            			if( queryEnd > sorghumIntervalEnd ) {
						            				queryEnd = sorghumIntervalEnd;
						            			}
//						            			System.out.println(samRecord.getReadName() +  " " + samRecord.getCigarString() + " " + samRecord.getAlignmentStart() + " " + samRecord.getAlignmentEnd() + " " + samRecord.getReadNegativeStrandFlag());
//							            		System.out.println("" + samRecord.getReadString().length() + "\t" + queryStart + " " + queryEnd);
							            		
								            	AlignedSorghumFragment alignedSorghumFragment = new AlignedSorghumFragment(samRecord.getReadName(), queryStart, queryEnd, strand);
								            	if( 0 == i ) {
								            		alignedSorghumFragments1.add(alignedSorghumFragment);
								            	}else {
								            		alignedSorghumFragments2.add(alignedSorghumFragment);
								            	}
						            		}
						            	}
						            }
								}
//					            System.out.println();
//					            System.out.println();
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
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum1.get(sorghumIntervalChr).charAt(position-1) != 'n' )  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum1.get(sorghumIntervalChr).charAt(position-1) != 'n')  {
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
					
					if( sorghumGenes.get(sorghumGene).getNumberOfCds() > 1 ) { // inner
						String sorghumIntervalChr = sorghumGenes.get(sorghumGene).getChromeSomeName();
						int sorghumIntervalStart = sorghumGenes.get(sorghumGene).getStart();
						int sorghumIntervalEnd = sorghumGenes.get(sorghumGene).getEnd();
						
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment>();
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment>();
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
							            			if( queryStart < sorghumIntervalStart ) {
							            				queryStart = sorghumIntervalStart;
							            			}
							            			if( queryEnd > sorghumIntervalEnd ) {
							            				queryEnd = sorghumIntervalEnd;
							            			}
							            			//System.out.println("" + samRecord.getReadString().length() + "\t" + queryStart + " " + queryEnd);
									            	AlignedSorghumFragment alignedSorghumFragment = new AlignedSorghumFragment(samRecord.getReadName(), queryStart, queryEnd, strand);
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
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n')  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum2.get(sorghumIntervalChr).charAt(position-1) != 'n')  {
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
						
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments1 = new ArrayList<AlignedSorghumFragment>();
						ArrayList<AlignedSorghumFragment> alignedSorghumFragments2 = new ArrayList<AlignedSorghumFragment>();
						for (int i=0; i<2; ++i) {
							int wantedStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() != sorghumGenes.get(sorghumGene).getStrand() ) {
								wantedStrand = 0;
							}
							
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+1;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeExtendLength;
							
							if ( 0 == wantedStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeExtendLength;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-1;
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
						            			if( queryStart < sorghumIntervalStart ) {
						            				queryStart = sorghumIntervalStart;
						            			}
						            			if( queryEnd > sorghumIntervalEnd ) {
						            				queryEnd = sorghumIntervalEnd;
						            			}
								            	AlignedSorghumFragment alignedSorghumFragment = new AlignedSorghumFragment(samRecord.getReadName(), queryStart, queryEnd, strand);
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
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments1 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum1.get(sorghumIntervalChr).charAt(position-1) != 'n')  {
									allPositions.add(position);
									positions1.add(position);
								}
							}
						}
						for( AlignedSorghumFragment alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position = alignedSorghumFragment.getStart(); position<=alignedSorghumFragment.getEnd(); ++position ) {
								if( sorghum1.get(sorghumIntervalChr).charAt(position-1) != 'n')  {
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
}

class AlignedSorghumFragment{
	private String chr;
	private int start;
	private int end;
	private int strand; //1 positive, 0 negative
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public int getStrand() {
		return strand;
	}
	public void setStrand(int strand) {
		this.strand = strand;
	}
	public AlignedSorghumFragment(String chr, int start, int end, int strand) {
		super();
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.strand = strand;
	}
}
