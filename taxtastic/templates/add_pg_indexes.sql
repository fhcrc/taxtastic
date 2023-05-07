{% set schema = schema or 'public' %}

CREATE INDEX ix_ranks_rank ON {{ schema }}.ranks (rank);

CREATE INDEX ix_nodes_tax_id ON {{ schema }}.nodes (tax_id);
CREATE INDEX ix_nodes_parent_id ON {{ schema }}.nodes (parent_id);
CREATE INDEX ix_nodes_source_id ON {{ schema }}.nodes (source_id);
CREATE INDEX ix_nodes_rank ON {{ schema }}.nodes (rank);

CREATE INDEX ix_names_is_primary ON {{ schema }}.names (is_primary);
CREATE INDEX ix_names_tax_id ON {{ schema }}.names (tax_id);
CREATE INDEX ix_names_tax_id_is_primary ON {{ schema }}.names (tax_id, is_primary);

CREATE INDEX ix_merged_old_tax_id ON {{ schema }}.merged (old_tax_id);
