
package us.kbase.genomefileutil;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: GenbankToGenomeParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "file",
    "genome_name",
    "workspace_name",
    "source",
    "taxon_wsname"
})
public class GenbankToGenomeParams {

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("file")
    private File file;
    @JsonProperty("genome_name")
    private String genomeName;
    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("source")
    private String source;
    @JsonProperty("taxon_wsname")
    private String taxonWsname;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("file")
    public File getFile() {
        return file;
    }

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("file")
    public void setFile(File file) {
        this.file = file;
    }

    public GenbankToGenomeParams withFile(File file) {
        this.file = file;
        return this;
    }

    @JsonProperty("genome_name")
    public String getGenomeName() {
        return genomeName;
    }

    @JsonProperty("genome_name")
    public void setGenomeName(String genomeName) {
        this.genomeName = genomeName;
    }

    public GenbankToGenomeParams withGenomeName(String genomeName) {
        this.genomeName = genomeName;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public GenbankToGenomeParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("source")
    public String getSource() {
        return source;
    }

    @JsonProperty("source")
    public void setSource(String source) {
        this.source = source;
    }

    public GenbankToGenomeParams withSource(String source) {
        this.source = source;
        return this;
    }

    @JsonProperty("taxon_wsname")
    public String getTaxonWsname() {
        return taxonWsname;
    }

    @JsonProperty("taxon_wsname")
    public void setTaxonWsname(String taxonWsname) {
        this.taxonWsname = taxonWsname;
    }

    public GenbankToGenomeParams withTaxonWsname(String taxonWsname) {
        this.taxonWsname = taxonWsname;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((("GenbankToGenomeParams"+" [file=")+ file)+", genomeName=")+ genomeName)+", workspaceName=")+ workspaceName)+", source=")+ source)+", taxonWsname=")+ taxonWsname)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
