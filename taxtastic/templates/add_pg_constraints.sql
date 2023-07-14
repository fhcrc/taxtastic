-- statements were created by pg_dump in this order
{% set schema = schema or 'public' %}

ALTER TABLE ONLY {{ schema }}.merged ADD CONSTRAINT merged_pkey PRIMARY KEY (old_tax_id);
ALTER TABLE ONLY {{ schema }}.names ADD CONSTRAINT names_pkey PRIMARY KEY (tax_id, tax_name, name_class);
ALTER TABLE ONLY {{ schema }}.nodes ADD CONSTRAINT nodes_pkey PRIMARY KEY (tax_id);
ALTER TABLE ONLY {{ schema }}.ranks ADD CONSTRAINT ranks_height_key UNIQUE (height);
ALTER TABLE ONLY {{ schema }}.ranks ADD CONSTRAINT ranks_pkey PRIMARY KEY (rank);
ALTER TABLE ONLY {{ schema }}.source ADD CONSTRAINT source_name_key UNIQUE (name);
ALTER TABLE ONLY {{ schema }}.source ADD CONSTRAINT source_pkey PRIMARY KEY (id);
ALTER TABLE ONLY {{ schema }}.merged ADD CONSTRAINT merged_new_tax_id_fkey FOREIGN KEY (new_tax_id) REFERENCES {{ schema }}.nodes(tax_id) ON DELETE CASCADE;
ALTER TABLE ONLY {{ schema }}.names ADD CONSTRAINT names_source_id_fkey FOREIGN KEY (source_id) REFERENCES {{ schema }}.source(id);
ALTER TABLE ONLY {{ schema }}.names ADD CONSTRAINT names_tax_id_fkey FOREIGN KEY (tax_id) REFERENCES {{ schema }}.nodes(tax_id) ON DELETE CASCADE;
ALTER TABLE ONLY {{ schema }}.nodes ADD CONSTRAINT nodes_rank_fkey FOREIGN KEY (rank) REFERENCES {{ schema }}.ranks(rank);
ALTER TABLE ONLY {{ schema }}.nodes ADD CONSTRAINT nodes_source_id_fkey FOREIGN KEY (source_id) REFERENCES {{ schema }}.source(id);
