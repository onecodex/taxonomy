//! Code related to handling of taxonomic ranks
use serde::{Deserialize, Serialize};
use std::str::FromStr;

use crate::{Result, TaxonomyError};

/// A taxonomic rank. For example, a species or phylum.
///
/// We use this instead of a String/&str to allow stricter type-checking
/// by forcing all taxonomic ranks to fall within the below categories
/// (this includes all current NCBI ranks and a few others, mostly ones
/// specific to zoology and botany).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TaxRank {
    Domain,
    Subdomain,
    Realm,
    Subrealm,
    Hyperkingdom,
    Superkingdom,
    Kingdom,
    Subkingdom,
    Infrakingdom,
    Parvkingdom,
    Superphylum,
    Phylum,
    Subphylum,
    Infraphylum,
    Microphylum,
    Superclass,
    Class,
    Subclass,
    Infraclass,
    Parvclass,
    Superdivision,
    Division,
    Subdivision,
    Infradivision,
    Superlegion,
    Legion,
    Sublegion,
    Infralegion,
    Supercohort,
    Cohort,
    Subcohort,
    Infracohort,
    Superorder,
    Gigaorder,
    Magnorder,
    Grandorder,
    Mirorder,
    SeriesFish,
    Order,
    Nanorder,
    Hypoorder,
    Minorder,
    Suborder,
    Infraorder,
    Parvorder,
    Section,
    Subsection,
    Gigafamily,
    Megafamily,
    Grandfamily,
    Hyperfamily,
    Superfamily,
    Epifamily,
    SeriesLepidoptera,
    GroupLepidoptera,
    Family,
    Subfamily,
    Infrafamily,
    Supertribe,
    Tribe,
    Subtribe,
    Infratribe,
    Genus,
    Subgenus,
    Series,
    SubseriesBotany,
    SpeciesGroup,
    SpeciesSubgroup,
    Species,
    Subspecies,
    Varietas,
    Subvarietas,
    Forma,
    Subforma,
    Cultivar,
    Breed,
    Strain,
    Individual,
    Clade,
    SeroGroup,
    Biotype,
    FormaSpecialis,
    Isolate,
    Serotype,
    Genotype,
    Morph,
    Pathogroup,
    // TODO: Unspecified prevents an auto-impl of Ord because it has no defined
    // place in the ordering (like a NaN) so we should manually derive out a
    // PartialOrd impl for TaxRank
    Unspecified,
}

impl TaxRank {
    /// Coverts a TaxRank into a one of the rank strings NCBI uses.
    /// Note that this doesn't handle ranks that are not used by the NCBI taxonomy.
    pub fn to_ncbi_rank(self) -> &'static str {
        match self {
            TaxRank::Superkingdom => "superkingdom",
            TaxRank::Kingdom => "kingdom",
            TaxRank::Subkingdom => "subkingdom",
            TaxRank::Superphylum => "superphylum",
            TaxRank::Phylum => "phylum",
            TaxRank::Subphylum => "subphylum",
            TaxRank::Superclass => "superclass",
            TaxRank::Class => "class",
            TaxRank::Subclass => "subclass",
            TaxRank::Infraclass => "infraclass",
            TaxRank::Cohort => "cohort",
            TaxRank::Subcohort => "subcohort",
            TaxRank::Superorder => "superorder",
            TaxRank::Order => "order",
            TaxRank::Suborder => "suborder",
            TaxRank::Infraorder => "infraorder",
            TaxRank::Parvorder => "parvorder",
            TaxRank::Superfamily => "superfamily",
            TaxRank::Family => "family",
            TaxRank::Subfamily => "subfamily",
            TaxRank::Tribe => "tribe",
            TaxRank::Subtribe => "subtribe",
            TaxRank::Genus => "genus",
            TaxRank::Series => "series",
            TaxRank::Subgenus => "subgenus",
            TaxRank::SpeciesGroup => "species group",
            TaxRank::SpeciesSubgroup => "species subgroup",
            TaxRank::Species => "species",
            TaxRank::Subspecies => "subspecies",
            TaxRank::Strain => "strain",
            TaxRank::Varietas => "varietas",
            TaxRank::Forma => "forma",
            TaxRank::Unspecified => "no rank",
            TaxRank::Clade => "clade",
            TaxRank::SeroGroup => "serogroup",
            TaxRank::Biotype => "biotype",
            TaxRank::FormaSpecialis => "forma specialis",
            TaxRank::Isolate => "isolate",
            TaxRank::Serotype => "serotype",
            TaxRank::Genotype => "genotype",
            TaxRank::Morph => "morph",
            TaxRank::Pathogroup => "pathogroup",
            // TODO: not sure if we want to manually coerce everything like this?
            _ => "no rank",
        }
    }
}

impl FromStr for TaxRank {
    type Err = TaxonomyError;

    fn from_str(s: &str) -> Result<Self> {
        // many of these synonyms (and the ranks themselves) were pulled from:
        // https://en.wikipedia.org/wiki/Taxonomic_rank
        match s.trim().to_lowercase().as_ref() {
            "domain" | "regio" => Ok(TaxRank::Domain),
            "subdomain" => Ok(TaxRank::Subdomain),
            "realm" => Ok(TaxRank::Realm),
            "subrealm" => Ok(TaxRank::Subrealm),
            "hyperkingdom" | "hyperregnum" => Ok(TaxRank::Hyperkingdom),
            "superkingdom" | "superregnum" => Ok(TaxRank::Superkingdom),
            "kingdom" | "regnum" => Ok(TaxRank::Kingdom),
            "subkingdom" | "subregnum" => Ok(TaxRank::Subkingdom),
            "infrakingdom" | "infraregnum" => Ok(TaxRank::Infrakingdom),
            "parvkingdom" | "parvregnum" => Ok(TaxRank::Parvkingdom),
            "superphylum" | "superphyla" => Ok(TaxRank::Superphylum),
            "phylum" | "phyla" => Ok(TaxRank::Phylum),
            "subphylum" | "subphyla" => Ok(TaxRank::Subphylum),
            "infraphylum" | "infraphyla" => Ok(TaxRank::Infraphylum),
            "microphylum" | "microphyla" => Ok(TaxRank::Microphylum),
            "superclass" => Ok(TaxRank::Superclass),
            "class" | "classis" => Ok(TaxRank::Class),
            "subclass" | "subclassis" => Ok(TaxRank::Subclass),
            "infraclass" => Ok(TaxRank::Infraclass),
            "parvclass" => Ok(TaxRank::Parvclass),
            // note "division" may be either a TaxRank::Division (zoological)
            // or synonum for TaxRank::Phylum (botanical) so we don't impl here
            "superlegion" => Ok(TaxRank::Superlegion),
            "legion" => Ok(TaxRank::Legion),
            "sublegion" => Ok(TaxRank::Sublegion),
            "infralegion" => Ok(TaxRank::Infralegion),
            "supercohort" => Ok(TaxRank::Supercohort),
            "cohort" => Ok(TaxRank::Cohort),
            "subcohort" => Ok(TaxRank::Subcohort),
            "infracohort" => Ok(TaxRank::Infracohort),
            "superorder" => Ok(TaxRank::Superorder),
            "gigaorder" => Ok(TaxRank::Gigaorder),
            "magnorder" => Ok(TaxRank::Magnorder),
            "grandorder" => Ok(TaxRank::Grandorder),
            "mirorder" => Ok(TaxRank::Mirorder),
            "order" | "ordo" => Ok(TaxRank::Order),
            "nanorder" => Ok(TaxRank::Nanorder),
            "hypoorder" => Ok(TaxRank::Hypoorder),
            "minorder" => Ok(TaxRank::Minorder),
            "suborder" | "subordo" => Ok(TaxRank::Suborder),
            "infraorder" => Ok(TaxRank::Infraorder),
            "parvorder" => Ok(TaxRank::Parvorder),
            "section" | "sectio" => Ok(TaxRank::Section),
            "subsection" => Ok(TaxRank::Subsection),
            "gigafamily" => Ok(TaxRank::Gigafamily),
            "megafamily" => Ok(TaxRank::Megafamily),
            "grandfamily" => Ok(TaxRank::Grandfamily),
            "hyperfamily" => Ok(TaxRank::Hyperfamily),
            "superfamily" => Ok(TaxRank::Superfamily),
            "epifamily" => Ok(TaxRank::Epifamily),
            "family" | "familia" => Ok(TaxRank::Family),
            "subfamily" => Ok(TaxRank::Subfamily),
            "infrafamily" => Ok(TaxRank::Infrafamily),
            "supertribe" => Ok(TaxRank::Supertribe),
            "tribe" | "subtribus" => Ok(TaxRank::Tribe),
            "subtribe" => Ok(TaxRank::Subtribe),
            "infratribe" => Ok(TaxRank::Infratribe),
            "genus" | "genera" => Ok(TaxRank::Genus),
            "subgenus" => Ok(TaxRank::Subgenus),
            "series" => Ok(TaxRank::Series),
            "species group" => Ok(TaxRank::SpeciesGroup),
            "species subgroup" => Ok(TaxRank::SpeciesSubgroup),
            "species" => Ok(TaxRank::Species),
            "subspecies" => Ok(TaxRank::Subspecies),
            "variety" | "varietas" => Ok(TaxRank::Varietas),
            "subvariety" | "subvarietas" => Ok(TaxRank::Subvarietas),
            "form" | "forma" => Ok(TaxRank::Forma),
            "subform" | "subforma" => Ok(TaxRank::Subforma),
            "cultivar" => Ok(TaxRank::Cultivar),
            "breed" => Ok(TaxRank::Breed),
            "strain" => Ok(TaxRank::Strain),
            "serogroup" => Ok(TaxRank::SeroGroup),
            "no rank" => Ok(TaxRank::Unspecified),
            "biotype" => Ok(TaxRank::Biotype),
            "clade" => Ok(TaxRank::Clade),
            "forma specialis" => Ok(TaxRank::FormaSpecialis),
            "isolate" => Ok(TaxRank::Isolate),
            "serotype" => Ok(TaxRank::Serotype),
            "genotype" => Ok(TaxRank::Genotype),
            "morph" => Ok(TaxRank::Morph),
            "pathogroup" => Ok(TaxRank::Pathogroup),
            _ => Err(TaxonomyError::UnrecognizedRank {
                rank: s.to_string(),
            }),
        }
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use crate::Result;

    use super::TaxRank;
    use super::TaxRank::*;

    static RANKS: &[super::TaxRank] = &[
        Domain,
        Subdomain,
        Realm,
        Subrealm,
        Hyperkingdom,
        Superkingdom,
        Kingdom,
        Subkingdom,
        Infrakingdom,
        Parvkingdom,
        Superphylum,
        Phylum,
        Subphylum,
        Infraphylum,
        Microphylum,
        Superclass,
        Class,
        Subclass,
        Infraclass,
        Parvclass,
        Superdivision,
        Division,
        Subdivision,
        Infradivision,
        Superlegion,
        Legion,
        Sublegion,
        Infralegion,
        Supercohort,
        Cohort,
        Subcohort,
        Infracohort,
        Superorder,
        Gigaorder,
        Magnorder,
        Grandorder,
        Mirorder,
        SeriesFish,
        Order,
        Nanorder,
        Hypoorder,
        Minorder,
        Suborder,
        Infraorder,
        Parvorder,
        Section,
        Subsection,
        Gigafamily,
        Megafamily,
        Grandfamily,
        Hyperfamily,
        Superfamily,
        Epifamily,
        SeriesLepidoptera,
        GroupLepidoptera,
        Family,
        Subfamily,
        Infrafamily,
        Supertribe,
        Tribe,
        Subtribe,
        Infratribe,
        Genus,
        Subgenus,
        Series,
        SubseriesBotany,
        SpeciesGroup,
        SpeciesSubgroup,
        Species,
        Subspecies,
        Varietas,
        Subvarietas,
        Forma,
        Subforma,
        Cultivar,
        Breed,
        Strain,
        Individual,
        Unspecified,
    ];

    #[test]
    fn test_ranks() {
        for rank in RANKS.iter() {
            let _ = rank.to_ncbi_rank();
        }
    }

    #[test]
    fn test_str_to_rank() -> Result<()> {
        for rank in RANKS.iter() {
            let _ = TaxRank::from_str(rank.to_ncbi_rank())?;
        }
        assert!(TaxRank::from_str("fake_data").is_err());
        Ok(())
    }
}
