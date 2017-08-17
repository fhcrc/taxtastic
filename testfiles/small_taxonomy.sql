-- Staphylococcaceae 90964
-- Gemella 1378

-- Keep only tax_ids in the lineage above these tax_ids and all of
-- their descendants.

.echo on

CREATE TEMPORARY TABLE temptab (tax_id text);

INSERT INTO temptab(tax_id)
WITH RECURSIVE a AS (
 SELECT tax_id, parent_id, rank
  FROM nodes
  WHERE tax_id in ('1378', '90964')
UNION ALL
 SELECT p.tax_id, p.parent_id, p.rank
  FROM a JOIN nodes p ON a.parent_id = p.tax_id
)
SELECT tax_id FROM a;

INSERT INTO temptab(tax_id)
WITH RECURSIVE descendants AS (
 SELECT tax_id, parent_id, rank
 FROM nodes
 WHERE tax_id in ('1378', '90964')
 UNION
 SELECT n.tax_id, n.parent_id, n.rank
 FROM nodes n
 INNER JOIN descendants d ON d.tax_id = n.parent_id
)
SELECT tax_id from descendants;

-- select names.tax_id, names.tax_name, nodes.*
-- from names
-- join nodes using(tax_id)
-- join ranks using(rank)
-- where tax_id in
-- (select tax_id from temptab)
-- and is_primary
-- order by height desc, tax_name;

delete from nodes
where tax_id not in
(select tax_id from temptab);

delete from names
where tax_id not in
(select tax_id from nodes);

delete from merged
where new_tax_id not in
(select tax_id from nodes);

vacuum;
