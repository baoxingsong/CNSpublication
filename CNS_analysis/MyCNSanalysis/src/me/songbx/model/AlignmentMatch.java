package me.songbx.model;


public class AlignmentMatch {
	private String species;
	private String chr;
	private int position;
	public String getSpecies() {
		return species;
	}
	public void setSpecies(String species) {
		this.species = species;
	}
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getPosition() {
		return position;
	}
	public void setPosition(int position) {
		this.position = position;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + position;
		result = prime * result + ((species == null) ? 0 : species.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AlignmentMatch other = (AlignmentMatch) obj;
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (position != other.position)
			return false;
		if (species == null) {
			if (other.species != null)
				return false;
		} else if (!species.equals(other.species))
			return false;
		return true;
	}
	public AlignmentMatch(String species, String chr, int position) {
		super();
		this.species = species;
		this.chr = chr;
		this.position = position;
	}
}
