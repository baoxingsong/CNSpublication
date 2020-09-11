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
import me.songbx.model.V3SnpsIdToV4Position;
import me.songbx.service.IfIntron;

public class CNSCDSOtherMAF {
	class SNP{
		private int position;
		private double maf;
		public int getPosition() {
			return position;
		}
		public void setPosition(int position) {
			this.position = position;
		}
		public double getMaf() {
			return maf;
		}
		public void setMaf(double maf) {
			this.maf = maf;
		}
		public SNP(int position, double maf) {
			this.position = position;
			this.maf = maf;
		}
	}
	
	public ArrayList<SNP> readMafFile( String mafFile, HashMap<String, Integer> snpIdToPosition ){
		ArrayList<SNP> snps = new ArrayList<SNP>();
		try {
        	File file = new File(mafFile);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
			
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = (tempString.trim()).split("\\s+");
				if( snpIdToPosition.containsKey( arrOfStr[1]) && arrOfStr[4].compareTo("NA")!=0  ) {
					SNP snp = new SNP(snpIdToPosition.get(arrOfStr[1]), Double.parseDouble(arrOfStr[4]));
					snps.add(snp);
				}else {
//					System.out.println("did not find " + arrOfStr[1]);
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		return snps;
	}
	
	public void readWigFile(String wigFile, HashMap<String, HashSet<Integer>> wigRecords) {
		try {
			File file = new File( wigFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String chr = "";
			Pattern p =Pattern.compile("chrom=(\\S+)");
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) == 't' ) {
					
				}else if( tempString.charAt(0) == 'v' ) {
					Matcher m = p.matcher(tempString);
					if( m.find() ) {
						chr = m.group(1);
					}
				}else {
					String[] arrOfStr = tempString.split("\\s+");
					if( wigRecords.containsKey(chr) ) {
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							wigRecords.get(chr).add(Integer.parseInt(arrOfStr[0]));
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public CNSCDSOtherMAF( String wigFile, String outputFile ){
		
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/Zea_mays.AGPv4.dna.toplevel_callable.fa").getChromoSomeHashMap();
		
		System.out.println("v3 to v4 reading begin");
		HashMap<String, Integer> snpIdToPosition = new V3SnpsIdToV4Position().getSnpIdToPosition();
		System.out.println("v3 to v4 reading done");
		HashMap<String, ArrayList<SNP>> snpHashMap = new HashMap<String, ArrayList<SNP>>();
		{//setup maf data structure begin
			System.out.println("allele frequency reading begin");
			ArrayList<SNP> snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr1.frq", snpIdToPosition);
			snpHashMap.put("1", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr2.frq", snpIdToPosition);
			snpHashMap.put("2", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr3.frq", snpIdToPosition);
			snpHashMap.put("3", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr4.frq", snpIdToPosition);
			snpHashMap.put("4", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr5.frq", snpIdToPosition);
			snpHashMap.put("5", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr6.frq", snpIdToPosition);
			snpHashMap.put("6", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr7.frq", snpIdToPosition);
			snpHashMap.put("7", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr8.frq", snpIdToPosition);
			snpHashMap.put("8", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr9.frq", snpIdToPosition);
			snpHashMap.put("9", snps);
			snps.clear();
			snps = readMafFile("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/chr10.frq", snpIdToPosition);
			snpHashMap.put("10", snps);
		}//setup maf data structure end
		System.out.println("allele frequency reading done");
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		System.out.println("gff file reading done");
		HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
		{
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			System.out.println("wig file reading begin");
			readWigFile(wigFile, cnsRecords);
		}
		System.out.println("wig file reading done");
		HashSet<Character> validatedDnaChars = new HashSet<Character>();
		validatedDnaChars.add('A');
		validatedDnaChars.add('T');
		validatedDnaChars.add('C');
		validatedDnaChars.add('G');
		validatedDnaChars.add('a');
		validatedDnaChars.add('t');
		validatedDnaChars.add('c');
		validatedDnaChars.add('g');
		try {
			PrintWriter outPut = new PrintWriter(outputFile);
    		
			for( String chr: snpHashMap.keySet() ) {
				for( SNP snp : snpHashMap.get(chr) ) {
					if ( validatedDnaChars.contains(maizeGenome.get(chr).charAt(snp.getPosition()-1)) ) {
						
						int elements = ifIntron.getElement(chr, snp.getPosition());
						if( elements == 2) {
							outPut.write(chr + "\t" + snp.getPosition() + "\t" + snp.getMaf() + "\tCDS\t0\n");
						}else {
							if (cnsRecords.containsKey(chr) && cnsRecords.get(chr).contains(snp.getPosition()) ){
								outPut.write(chr + "\t" + snp.getPosition() + "\t" + snp.getMaf() + "\tCNS\t0\n");
							} else if(0==elements) {
								outPut.write(chr + "\t" + snp.getPosition() + "\t" + snp.getMaf() + "\tintergenetic\t" + ifIntron.getMinidistanceWithCds(chr, snp.getPosition()).getMiniDistance() + "\n");
							} else if(3==elements) {
								outPut.write(chr + "\t" + snp.getPosition() + "\t" + snp.getMaf() + "\tUTR\t" + ifIntron.getMinidistanceWithCds(chr, snp.getPosition()).getMiniDistance() + "\n");
							} else {
								outPut.write(chr + "\t" + snp.getPosition() + "\t" + snp.getMaf() + "\tintron\t" + ifIntron.getMinidistanceWithCds(chr, snp.getPosition()).getMiniDistance() + "\n");
							}
						}
					}
				}
			}
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		new CNSCDSOtherMAF("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.wig", "/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/pan_and_cns_depth.maf");
		new CNSCDSOtherMAF("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_depth.wig", "/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/core_and_cns_depth.maf");
		new CNSCDSOtherMAF("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_mpileup.wig", "/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/pan_and_cns_mpileup.maf");
		new CNSCDSOtherMAF("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_mpileup.wig", "/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/core_and_cns_mpileup.maf");
//		new CNSCDSOtherMAF(args[0], args[1]);
		// java -jar -Xmx25g CNSCDSOtherMAF.jar /media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.bw /media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/pan_and_cns_depth.maf
		// java -jar -Xmx25g CNSCDSOtherMAF.jar /media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_depth.bw /media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/core_and_cns_depth.maf
		// java -jar -Xmx25g CNSCDSOtherMAF.jar /media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_mpileup.bw /media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/pan_and_cns_mpileup.maf
		// java -jar -Xmx25g CNSCDSOtherMAF.jar /media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_mpileup.bw /media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/coreAndCns/core_and_cns_mpileup.maf
	}
}
