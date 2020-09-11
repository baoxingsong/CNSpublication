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
import me.songbx.model.NearByGene;
import me.songbx.service.IfIntron;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class GenesWithCNSDetected {
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
	
	public GenesWithCNSDetected( String wigFile, String outputFile ){
		
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/Zea_mays.AGPv4.dna.toplevel_callable.fa").getChromoSomeHashMap();
		
		
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
			
			ExecutorService myExecutor = Executors.newFixedThreadPool(10);
			NearyGene genesInRange1 = new NearyGene();
			NearyGene genesInRange2 = new NearyGene();
			for( String chr: cnsRecords.keySet() ) {
				for( Integer position : cnsRecords.get(chr) ) {
					if ( validatedDnaChars.contains(maizeGenome.get(chr).charAt(position-1)) ) {
						myExecutor.execute(new RunItParallele( ifIntron, chr, position, genesInRange1, genesInRange2));
					}
				}
			}
			
			myExecutor.shutdown();
			try {
				myExecutor.awaitTermination(2000, TimeUnit.HOURS);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			
			PrintWriter outPut = new PrintWriter(outputFile);
			for( String gene: genesInRange1.getGenesInRange().keySet() ) {
				outPut.write(gene + "\t" + genesInRange1.getGenesInRange().get(gene) + "\n");
			}
			
			for( String gene: genesInRange2.getGenesInRange().keySet() ) {
				outPut.write(gene + "\t" + genesInRange2.getGenesInRange().get(gene) + "\n");
			}
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	class RunItParallele implements Runnable{
		
		private IfIntron ifIntron;
		private String chr;
		private int position;
		private NearyGene genesInRange1;
		private NearyGene genesInRange2;
		public RunItParallele(IfIntron ifIntron, String chr, int position, NearyGene genesInRange1, NearyGene genesInRange2) {
			this.ifIntron=ifIntron;
			this.chr=chr;
			this.position=position;
			this.genesInRange1=genesInRange1;
			this.genesInRange2=genesInRange2;
		}

		@Override
		public void run() {
			int elements = ifIntron.getElement(chr, position);
			if(3==elements || 0==elements) {
				NearByGene nearByGene = ifIntron.getMinidistanceWithCds(chr, position);
				if( nearByGene.getMiniDistance()>0 ) {
					if (  (!genesInRange1.containsKey(nearByGene.getGene())) || (genesInRange1.get(nearByGene.getGene()) > nearByGene.getMiniDistance())    )  {
						genesInRange1.put(nearByGene.getGene(), nearByGene.getMiniDistance());
					}
				} else {
					if (  (!genesInRange2.containsKey(nearByGene.getGene())) || (genesInRange2.get(nearByGene.getGene()) < nearByGene.getMiniDistance())    )  {
						genesInRange2.put(nearByGene.getGene(), nearByGene.getMiniDistance());
					}
				}
			}		
		}
	}
	
	class NearyGene{
		private HashMap<String, Integer> genesInRange = new HashMap<String, Integer>();
		public synchronized void put( String key, Integer value) {
			genesInRange.put(key, value);
		}
		public synchronized Integer get( String key) {
			return genesInRange.get(key);
		}
		public synchronized boolean containsKey( String key) {
			return genesInRange.containsKey(key);
		}
		public HashMap<String, Integer> getGenesInRange() {
			return genesInRange;
		}
		public void setGenesInRange(HashMap<String, Integer> genesInRange) {
			this.genesInRange = genesInRange;
		}
	}
	
	public static void main(String[] args) {
			new GenesWithCNSDetected("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.wig", "/media/bs674/1_8t/AndCns/genewithCNS/pan_and_cns_depth.maf");
			new GenesWithCNSDetected("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/core_and_cns_depth.wig", "/media/bs674/1_8t/AndCns/genewithCNS/core_and_cns_depth.maf");

	}
}
