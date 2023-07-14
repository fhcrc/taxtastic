{% set schema = schema or 'public' %}

DROP INDEX ix_ranks_rank;
DROP INDEX ix_nodes_tax_id;
DROP INDEX ix_nodes_parent_id;
DROP INDEX ix_nodes_source_id;
DROP INDEX ix_nodes_rank;
DROP INDEX ix_names_is_primary;
DROP INDEX ix_names_tax_id;
DROP INDEX ix_names_tax_id_is_primary;
DROP INDEX ix_merged_old_tax_id;

ALTER TABLE ONLY {{ schema }}.merged DROP CONSTRAINT merged_new_tax_id_fkey;
ALTER TABLE ONLY {{ schema }}.merged DROP CONSTRAINT merged_pkey;

ALTER TABLE ONLY {{ schema }}.names DROP CONSTRAINT names_pkey;
ALTER TABLE ONLY {{ schema }}.names DROP CONSTRAINT names_source_id_fkey;
ALTER TABLE ONLY {{ schema }}.names DROP CONSTRAINT names_tax_id_fkey;

ALTER TABLE ONLY {{ schema }}.nodes DROP CONSTRAINT nodes_pkey;
ALTER TABLE ONLY {{ schema }}.nodes DROP CONSTRAINT nodes_rank_fkey;
ALTER TABLE ONLY {{ schema }}.nodes DROP CONSTRAINT nodes_source_id_fkey;

ALTER TABLE ONLY {{ schema }}.ranks DROP CONSTRAINT ranks_height_key;
ALTER TABLE ONLY {{ schema }}.ranks DROP CONSTRAINT ranks_pkey;

ALTER TABLE ONLY {{ schema }}.source DROP CONSTRAINT source_name_key;
ALTER TABLE ONLY {{ schema }}.source DROP CONSTRAINT source_pkey;
