package me.songbx.impl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import me.songbx.model.V3SnpsIdToV4Position;

public class LeadingSnps {
	class LeadingeQtl{
		private String chr;
		private double p_value;
		private String snpId;
		public String getSnpId() {
			return snpId;
		}
		public void setSnpId(String snpId) {
			this.snpId = snpId;
		}
		public String getChr() {
			return chr;
		}
		public void setChr(String chr) {
			this.chr = chr;
		}
		public double getP_value() {
			return p_value;
		}
		public void setP_value(double p_value) {
			this.p_value = p_value;
		}
		public LeadingeQtl(String chr, double p_value, String snpId) {
			super();
			this.chr = chr;
			this.p_value = p_value;
			this.snpId = snpId;
		}
	}
	private HashMap<String, HashSet<Integer>> leadingSnps;
	public HashMap<String, HashSet<Integer>> getLeadingSnps() {
		return leadingSnps;
	}
	public void setLeadingSnps(HashMap<String, HashSet<Integer>> leadingSnps) {
		this.leadingSnps = leadingSnps;
	}
	private void readeQTl(ArrayList<String> eQtlsFiles ) {
		HashMap<String, LeadingeQtl> leadingeQtls = new HashMap<String, LeadingeQtl>();
		for( String eQtlsFile : eQtlsFiles ) {
			System.out.println("reading:" + eQtlsFile);
			try {
				File file = new File( eQtlsFile );
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); //skip the first line
				while ((tempString = reader.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+");
					double p =  Double.parseDouble(arrOfStr[6]);
					String snpid = arrOfStr[1];
					snpid = snpid.replace("S", "");
					snpid = snpid.replace("_", "-");
					if( (!leadingeQtls.containsKey(arrOfStr[0]))  ||  leadingeQtls.get(arrOfStr[0]).getP_value() > p) {
						LeadingeQtl leadingeQtl = new LeadingeQtl(arrOfStr[2], p, snpid);
						leadingeQtls.put(arrOfStr[0], leadingeQtl);
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		HashMap<String, Integer> snpIdToPosition = new V3SnpsIdToV4Position().getSnpIdToPosition();
		
		for( String gene : leadingeQtls.keySet() ) {
			//System.out.println("line 114 " + gene + "\t" + leadingeQtls.get(gene).getSnpId());
			if( !leadingSnps.containsKey(leadingeQtls.get(gene).getChr())  ) {
				leadingSnps.put(leadingeQtls.get(gene).getChr(), new HashSet<Integer>());
			}
			if( snpIdToPosition.containsKey(leadingeQtls.get(gene).getSnpId()) ) {
				leadingSnps.get(leadingeQtls.get(gene).getChr()).add(snpIdToPosition.get(leadingeQtls.get(gene).getSnpId()));
				//System.out.println("line 119 " + snpIdToPosition.get(leadingeQtls.get(gene).getSnpId()));
			}
		}
	}
	private void generate() {
		{ //ROOT
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GRoot%c10");
			readeQTl(eQtlsFiles);
		}
		{ //SHOOT
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/GShoot%c10");
			readeQTl(eQtlsFiles);
		}
		{ //Kern
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/Kern%c10");
			readeQTl(eQtlsFiles);
		}
		{ //L3Base
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Base%c10");
			readeQTl(eQtlsFiles);
		}
		{ //L3Tip
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/L3Tip%c10");
			readeQTl(eQtlsFiles);
		}
		{ //LMAD
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAD%c10");
			readeQTl(eQtlsFiles);
		}
		{ //LMAN
			ArrayList<String> eQtlsFiles = new ArrayList<String>();
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c1");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c2");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c3");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c4");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c5");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c6");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c7");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c8");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c9");
			eQtlsFiles.add("/media/bs674/pan_and_non_asse/eQTL/LMAN%c10");
			readeQTl(eQtlsFiles);
		}
	}
	public HashMap<String, HashSet<Integer>> readSummary() {
		try {
			File file = new File( "/media/bs674/pan_and_non_asse/eQTL/leading_snps" );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = reader.readLine(); //skip the first line
			while ((tempString = reader.readLine()) != null) {
//				System.out.println(tempString);
				String[] arrOfStr = tempString.split("\\s+");
				if( !leadingSnps.containsKey(arrOfStr[0])  ) {
					leadingSnps.put(arrOfStr[0], new HashSet<Integer>());
				}
				leadingSnps.get(arrOfStr[0]).add(Integer.parseInt(arrOfStr[1]));
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return leadingSnps;
	}
	public LeadingSnps(){
		leadingSnps = new HashMap<String, HashSet<Integer>>();
		
	}
	// read eQTL result end
	public static void main( String argv[] ) {
		LeadingSnps leadingSnps = new LeadingSnps();
		
		leadingSnps.generate();
		HashMap<String, HashSet<Integer>> snps = leadingSnps.getLeadingSnps();
		try {
			PrintWriter outPut = new PrintWriter("/media/bs674/pan_and_non_asse/eQTL/leading_snps");
			for( String chr : snps.keySet() ) {
				for ( int position : snps.get(chr) ) {
					outPut.println(chr+"\t"+position);
				}
			}
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		leadingSnps.readSummary();
	}
}
