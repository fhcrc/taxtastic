{% for table in ['ranks', 'source', 'nodes', 'names', 'merged'] %}
{% if dialect == 'postgresql' %}
drop table if exists "{{ table }}" cascade;
{% else %}
drop table if exists "{{ table }}";
{% endif %}
{% endfor %}

CREATE TABLE ranks (
        rank TEXT NOT NULL,
        height INTEGER NOT NULL--,
        -- PRIMARY KEY (rank),
        -- UNIQUE (height)
);

CREATE TABLE source (
        id {% if dialect == 'postgresql' %}SERIAL{% else %}INTEGER{% endif %} NOT NULL,
        name TEXT,
        description TEXT,
        PRIMARY KEY (id),
        UNIQUE (name)
);

CREATE TABLE nodes (
        tax_id TEXT NOT NULL,
        parent_id TEXT,
        rank TEXT,
        embl_code TEXT,
        division_id TEXT,
        source_id INTEGER,
        is_valid BOOLEAN--,
        -- PRIMARY KEY (tax_id),
        -- FOREIGN KEY(rank) REFERENCES ranks (rank),
        -- FOREIGN KEY(source_id) REFERENCES source (id)
);

-- CREATE INDEX ix_nodes_parent_id ON nodes (parent_id);
CREATE TABLE names (
        tax_id TEXT NOT NULL,
        tax_name TEXT NOT NULL,
        unique_name TEXT,
        name_class TEXT NOT NULL,
        source_id INTEGER,
        is_primary BOOLEAN,
        is_classified BOOLEAN--,
        -- PRIMARY KEY (tax_id, tax_name, name_class),
        -- FOREIGN KEY(tax_id) REFERENCES nodes (tax_id) ON DELETE CASCADE,
        -- FOREIGN KEY(source_id) REFERENCES source (id)
);

-- CREATE INDEX ix_names_tax_id_is_primary ON names (tax_id, is_primary);

CREATE TABLE merged (
        old_tax_id TEXT NOT NULL,
        new_tax_id TEXT--,
        -- PRIMARY KEY (old_tax_id),
        -- FOREIGN KEY(new_tax_id) REFERENCES nodes (tax_id) ON DELETE CASCADE
);
-- CREATE INDEX ix_merged_old_tax_id ON merged (old_tax_id);

