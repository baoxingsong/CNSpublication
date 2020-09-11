package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import me.songbx.service.IfIntron;

public class IntronRegionMissingRatioAndMaf {

	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/intronMafAndMissing");
			IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
			System.out.println("gff reading done");
			for ( int i=1; i<=10; ++i ) {
				System.out.println("chr:" + i);
				HashMap<String, Integer> snpIdKeyPositionValue = new HashMap<String, Integer>();
				File file0 = new File("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr"+i+".vcf");
				BufferedReader reader0 = null;
				try {
					reader0 = new BufferedReader(new FileReader(file0));
		            String tempString = reader0.readLine(); // skip the header line
					while ((tempString = reader0.readLine()) != null) {
						if( tempString.charAt(0) != '#' ) {
							String[] currencies = tempString.split("\\s+");
							snpIdKeyPositionValue.put(currencies[2], Integer.parseInt(currencies[1]));
						}
			        }
					reader0.close();
					System.out.println("vcf file reading done " + snpIdKeyPositionValue.size());
		        }catch (IOException e) {
		            e.printStackTrace();
		        } finally {
		            if (reader0 != null) {
		                try {
		                	reader0.close();
		                } catch (IOException e1) {
		                	e1.getStackTrace();
		                }
		            }
		        }
				HashMap<String, String> missingRatio = new HashMap<String, String>();
				
				File file1 = new File("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/chr"+i+".lmiss");
				BufferedReader reader1 = null;
		        try {
		        	reader1 = new BufferedReader(new FileReader(file1));
		            String tempString = reader1.readLine(); // skip the header line
					while ((tempString = reader1.readLine()) != null) {
						tempString = tempString.trim();
						String[] currencies = tempString.split("\\s+");
						if ( currencies.length == 5 ) {
							if( 1 == ifIntron.getElement(Integer.toString(i), snpIdKeyPositionValue.get(currencies[1])) ) {
								missingRatio.put(currencies[1], currencies[4]);
							}
						}else {
							System.out.println(tempString + "\t\t\t\t" + currencies.length);
						}
			        }
					reader1.close();
		        }catch (IOException e) {
		            e.printStackTrace();
		        } finally {
		            if (reader1 != null) {
		                try {
		                    reader1.close();
		                } catch (IOException e1) {
		                	e1.getStackTrace();
		                }
		            }
		        }
		        System.out.println("msing ratio reading done " + missingRatio.size());
		        File file2 = new File("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/chr"+i+".frq");
				BufferedReader reader2 = null;
		        try {
		        	reader2 = new BufferedReader(new FileReader(file2));
		            String tempString = reader2.readLine(); // skip the header line
					while ((tempString = reader2.readLine()) != null) {
						tempString = tempString.trim();
						String[] currencies = tempString.split("\\s+");
						if ( currencies.length == 6 ) {
							if( 1 == ifIntron.getElement(Integer.toString(i), snpIdKeyPositionValue.get(currencies[1])) ) {
								outPut.println( currencies[1] + "\t" + Integer.toString(i) + "\t" + Integer.toString(snpIdKeyPositionValue.get(currencies[1])) + "\t" + missingRatio.get(currencies[1]) + "\t" + currencies[4] );
							}
						}
			        }
					reader2.close();
		        }catch (IOException e) {
		            e.printStackTrace();
		        } finally {
		            if (reader2 != null) {
		                try {
		                    reader2.close();
		                } catch (IOException e1) {
		                	e1.getStackTrace();
		                }
		            }
		        }
			}
			outPut.close();
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}finally {
            if (outPut != null) {
                outPut.close();
            }
        }
	}
}
