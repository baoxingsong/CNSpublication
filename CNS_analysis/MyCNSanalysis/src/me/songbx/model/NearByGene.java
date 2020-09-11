package me.songbx.model;

public class NearByGene {
	int miniDistance;
	String gene;
	public int getMiniDistance() {
		return miniDistance;
	}
	public void setMiniDistance(int miniDistance) {
		this.miniDistance = miniDistance;
	}
	public String getGene() {
		return gene;
	}
	public void setGene(String gene) {
		this.gene = gene;
	}
	public NearByGene(int miniDistance, String gene) {
		super();
		this.miniDistance = miniDistance;
		this.gene = gene;
	}
}
