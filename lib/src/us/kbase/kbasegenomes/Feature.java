
package us.kbase.kbasegenomes;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import us.kbase.common.service.Tuple2;
import us.kbase.common.service.Tuple4;


/**
 * <p>Original spec-file type: Feature</p>
 * <pre>
 * Structure for a single CDS encoding ???gene??? of a genome
 * ONLY PUT GENES THAT HAVE A CORRESPONDING CDS IN THIS ARRAY
 *     NOTE: Sequence is optional. Ideally we can keep it in here, but
 * Recognize due to space constraints another solution may be needed.
 *     We may want to add additional fields for other CDM functions
 *     (e.g., atomic regulons, coexpressed fids, co_occurring fids,...)
 *     protein_translation_length and protein_translation are
 *     for longest coded protein (representative protein for splice variants)
 *     NOTE: New Aliases field definitely breaks compatibility.
 *           As Does Function.
 *     flags are flag fields in GenBank format. This will be a controlled vocabulary.
 *     Initially Acceptable values are pseudo, ribosomal_slippage, and trans_splicing
 *     Md5 is the md5 of dna_sequence.
 *     @optional functions ontology_terms note protein_translation mrnas flags warnings
 *     @optional inference_data dna_sequence aliases db_xrefs children functional_descriptions
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "id",
    "location",
    "functions",
    "functional_descriptions",
    "ontology_terms",
    "note",
    "md5",
    "protein_translation",
    "protein_translation_length",
    "cdss",
    "mrnas",
    "children",
    "flags",
    "warnings",
    "inference_data",
    "dna_sequence",
    "dna_sequence_length",
    "aliases",
    "db_xrefs"
})
public class Feature {

    @JsonProperty("id")
    private java.lang.String id;
    @JsonProperty("location")
    private List<Tuple4 <String, Long, String, Long>> location;
    @JsonProperty("functions")
    private List<String> functions;
    @JsonProperty("functional_descriptions")
    private List<String> functionalDescriptions;
    @JsonProperty("ontology_terms")
    private Map<String, Map<String, List<Long>>> ontologyTerms;
    @JsonProperty("note")
    private java.lang.String note;
    @JsonProperty("md5")
    private java.lang.String md5;
    @JsonProperty("protein_translation")
    private java.lang.String proteinTranslation;
    @JsonProperty("protein_translation_length")
    private java.lang.Long proteinTranslationLength;
    @JsonProperty("cdss")
    private List<String> cdss;
    @JsonProperty("mrnas")
    private List<String> mrnas;
    @JsonProperty("children")
    private List<String> children;
    @JsonProperty("flags")
    private List<String> flags;
    @JsonProperty("warnings")
    private List<String> warnings;
    @JsonProperty("inference_data")
    private List<InferenceInfo> inferenceData;
    @JsonProperty("dna_sequence")
    private java.lang.String dnaSequence;
    @JsonProperty("dna_sequence_length")
    private java.lang.Long dnaSequenceLength;
    @JsonProperty("aliases")
    private List<Tuple2 <String, String>> aliases;
    @JsonProperty("db_xrefs")
    private List<Tuple2 <String, String>> dbXrefs;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("id")
    public java.lang.String getId() {
        return id;
    }

    @JsonProperty("id")
    public void setId(java.lang.String id) {
        this.id = id;
    }

    public Feature withId(java.lang.String id) {
        this.id = id;
        return this;
    }

    @JsonProperty("location")
    public List<Tuple4 <String, Long, String, Long>> getLocation() {
        return location;
    }

    @JsonProperty("location")
    public void setLocation(List<Tuple4 <String, Long, String, Long>> location) {
        this.location = location;
    }

    public Feature withLocation(List<Tuple4 <String, Long, String, Long>> location) {
        this.location = location;
        return this;
    }

    @JsonProperty("functions")
    public List<String> getFunctions() {
        return functions;
    }

    @JsonProperty("functions")
    public void setFunctions(List<String> functions) {
        this.functions = functions;
    }

    public Feature withFunctions(List<String> functions) {
        this.functions = functions;
        return this;
    }

    @JsonProperty("functional_descriptions")
    public List<String> getFunctionalDescriptions() {
        return functionalDescriptions;
    }

    @JsonProperty("functional_descriptions")
    public void setFunctionalDescriptions(List<String> functionalDescriptions) {
        this.functionalDescriptions = functionalDescriptions;
    }

    public Feature withFunctionalDescriptions(List<String> functionalDescriptions) {
        this.functionalDescriptions = functionalDescriptions;
        return this;
    }

    @JsonProperty("ontology_terms")
    public Map<String, Map<String, List<Long>>> getOntologyTerms() {
        return ontologyTerms;
    }

    @JsonProperty("ontology_terms")
    public void setOntologyTerms(Map<String, Map<String, List<Long>>> ontologyTerms) {
        this.ontologyTerms = ontologyTerms;
    }

    public Feature withOntologyTerms(Map<String, Map<String, List<Long>>> ontologyTerms) {
        this.ontologyTerms = ontologyTerms;
        return this;
    }

    @JsonProperty("note")
    public java.lang.String getNote() {
        return note;
    }

    @JsonProperty("note")
    public void setNote(java.lang.String note) {
        this.note = note;
    }

    public Feature withNote(java.lang.String note) {
        this.note = note;
        return this;
    }

    @JsonProperty("md5")
    public java.lang.String getMd5() {
        return md5;
    }

    @JsonProperty("md5")
    public void setMd5(java.lang.String md5) {
        this.md5 = md5;
    }

    public Feature withMd5(java.lang.String md5) {
        this.md5 = md5;
        return this;
    }

    @JsonProperty("protein_translation")
    public java.lang.String getProteinTranslation() {
        return proteinTranslation;
    }

    @JsonProperty("protein_translation")
    public void setProteinTranslation(java.lang.String proteinTranslation) {
        this.proteinTranslation = proteinTranslation;
    }

    public Feature withProteinTranslation(java.lang.String proteinTranslation) {
        this.proteinTranslation = proteinTranslation;
        return this;
    }

    @JsonProperty("protein_translation_length")
    public java.lang.Long getProteinTranslationLength() {
        return proteinTranslationLength;
    }

    @JsonProperty("protein_translation_length")
    public void setProteinTranslationLength(java.lang.Long proteinTranslationLength) {
        this.proteinTranslationLength = proteinTranslationLength;
    }

    public Feature withProteinTranslationLength(java.lang.Long proteinTranslationLength) {
        this.proteinTranslationLength = proteinTranslationLength;
        return this;
    }

    @JsonProperty("cdss")
    public List<String> getCdss() {
        return cdss;
    }

    @JsonProperty("cdss")
    public void setCdss(List<String> cdss) {
        this.cdss = cdss;
    }

    public Feature withCdss(List<String> cdss) {
        this.cdss = cdss;
        return this;
    }

    @JsonProperty("mrnas")
    public List<String> getMrnas() {
        return mrnas;
    }

    @JsonProperty("mrnas")
    public void setMrnas(List<String> mrnas) {
        this.mrnas = mrnas;
    }

    public Feature withMrnas(List<String> mrnas) {
        this.mrnas = mrnas;
        return this;
    }

    @JsonProperty("children")
    public List<String> getChildren() {
        return children;
    }

    @JsonProperty("children")
    public void setChildren(List<String> children) {
        this.children = children;
    }

    public Feature withChildren(List<String> children) {
        this.children = children;
        return this;
    }

    @JsonProperty("flags")
    public List<String> getFlags() {
        return flags;
    }

    @JsonProperty("flags")
    public void setFlags(List<String> flags) {
        this.flags = flags;
    }

    public Feature withFlags(List<String> flags) {
        this.flags = flags;
        return this;
    }

    @JsonProperty("warnings")
    public List<String> getWarnings() {
        return warnings;
    }

    @JsonProperty("warnings")
    public void setWarnings(List<String> warnings) {
        this.warnings = warnings;
    }

    public Feature withWarnings(List<String> warnings) {
        this.warnings = warnings;
        return this;
    }

    @JsonProperty("inference_data")
    public List<InferenceInfo> getInferenceData() {
        return inferenceData;
    }

    @JsonProperty("inference_data")
    public void setInferenceData(List<InferenceInfo> inferenceData) {
        this.inferenceData = inferenceData;
    }

    public Feature withInferenceData(List<InferenceInfo> inferenceData) {
        this.inferenceData = inferenceData;
        return this;
    }

    @JsonProperty("dna_sequence")
    public java.lang.String getDnaSequence() {
        return dnaSequence;
    }

    @JsonProperty("dna_sequence")
    public void setDnaSequence(java.lang.String dnaSequence) {
        this.dnaSequence = dnaSequence;
    }

    public Feature withDnaSequence(java.lang.String dnaSequence) {
        this.dnaSequence = dnaSequence;
        return this;
    }

    @JsonProperty("dna_sequence_length")
    public java.lang.Long getDnaSequenceLength() {
        return dnaSequenceLength;
    }

    @JsonProperty("dna_sequence_length")
    public void setDnaSequenceLength(java.lang.Long dnaSequenceLength) {
        this.dnaSequenceLength = dnaSequenceLength;
    }

    public Feature withDnaSequenceLength(java.lang.Long dnaSequenceLength) {
        this.dnaSequenceLength = dnaSequenceLength;
        return this;
    }

    @JsonProperty("aliases")
    public List<Tuple2 <String, String>> getAliases() {
        return aliases;
    }

    @JsonProperty("aliases")
    public void setAliases(List<Tuple2 <String, String>> aliases) {
        this.aliases = aliases;
    }

    public Feature withAliases(List<Tuple2 <String, String>> aliases) {
        this.aliases = aliases;
        return this;
    }

    @JsonProperty("db_xrefs")
    public List<Tuple2 <String, String>> getDbXrefs() {
        return dbXrefs;
    }

    @JsonProperty("db_xrefs")
    public void setDbXrefs(List<Tuple2 <String, String>> dbXrefs) {
        this.dbXrefs = dbXrefs;
    }

    public Feature withDbXrefs(List<Tuple2 <String, String>> dbXrefs) {
        this.dbXrefs = dbXrefs;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((((((((((((((((((("Feature"+" [id=")+ id)+", location=")+ location)+", functions=")+ functions)+", functionalDescriptions=")+ functionalDescriptions)+", ontologyTerms=")+ ontologyTerms)+", note=")+ note)+", md5=")+ md5)+", proteinTranslation=")+ proteinTranslation)+", proteinTranslationLength=")+ proteinTranslationLength)+", cdss=")+ cdss)+", mrnas=")+ mrnas)+", children=")+ children)+", flags=")+ flags)+", warnings=")+ warnings)+", inferenceData=")+ inferenceData)+", dnaSequence=")+ dnaSequence)+", dnaSequenceLength=")+ dnaSequenceLength)+", aliases=")+ aliases)+", dbXrefs=")+ dbXrefs)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
