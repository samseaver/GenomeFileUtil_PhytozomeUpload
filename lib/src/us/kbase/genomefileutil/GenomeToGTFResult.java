
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
 * <p>Original spec-file type: GenomeToGTFResult</p>
 * <pre>
 * from_cache is 1 if the file already exists and was just returned, 0 if
 * the file was generated during this call.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "gtf_file",
    "from_cache"
})
public class GenomeToGTFResult {

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gtf_file")
    private File gtfFile;
    @JsonProperty("from_cache")
    private Long fromCache;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gtf_file")
    public File getGtfFile() {
        return gtfFile;
    }

    /**
     * <p>Original spec-file type: File</p>
     * 
     * 
     */
    @JsonProperty("gtf_file")
    public void setGtfFile(File gtfFile) {
        this.gtfFile = gtfFile;
    }

    public GenomeToGTFResult withGtfFile(File gtfFile) {
        this.gtfFile = gtfFile;
        return this;
    }

    @JsonProperty("from_cache")
    public Long getFromCache() {
        return fromCache;
    }

    @JsonProperty("from_cache")
    public void setFromCache(Long fromCache) {
        this.fromCache = fromCache;
    }

    public GenomeToGTFResult withFromCache(Long fromCache) {
        this.fromCache = fromCache;
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
        return ((((((("GenomeToGTFResult"+" [gtfFile=")+ gtfFile)+", fromCache=")+ fromCache)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
