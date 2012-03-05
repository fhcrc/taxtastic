import sys; sys.path.append('../')
from taxtastic.lonely import *
import cStringIO
import os
from sqlalchemy import create_engine
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi


def test_trees():
    t = Tree(3, taxname="boris")(
          Tree(5, taxname="boris meep"),
          Tree(6, taxname="boris hilda"))
    assert t.key == 3
    assert t.taxname == 'boris'
    assert t.children[0].children == []
    assert t.children[0].taxname == "boris meep"
    assert t.children[0].parent == t
    assert t.descendents[5] == t.children[0]
    t2 = Tree(1)(
             Tree(3)(
                 Tree(4)(
                     Tree(5)(
                         Tree(7)),
                     Tree(6))))
    assert t2.descendents.keys() == [1,3,4,5,6,7]
    assert t2.children[0].descendents.keys() == [3,4,5,6,7]


taxtable = """"tax_id","parent_id","rank","tax_name","root","below_root","superkingdom","superphylum","phylum","class","subclass","order","below_order","suborder","family","genus","species"
"1","1","root","root","1","","","","","","","","","","","",""
"131567","1","below_root","cellular organisms","1","131567","","","","","","","","","","",""
"2","131567","superkingdom","Bacteria","1","131567","2","","","","","","","","","",""
"68336","2","superphylum","Bacteroidetes/Chlorobi group","1","131567","2","68336","","","","","","","","",""
"201174","2","phylum","Actinobacteria","1","131567","2","","201174","","","","","","","",""
"976","68336","phylum","Bacteroidetes","1","131567","2","68336","976","","","","","","","",""
"200930","2","phylum","Deferribacteres","1","131567","2","","200930","","","","","","","",""
"1239","2","phylum","Firmicutes","1","131567","2","","1239","","","","","","","",""
"32066","2","phylum","Fusobacteria","1","131567","2","","32066","","","","","","","",""
"1224","2","phylum","Proteobacteria","1","131567","2","","1224","","","","","","","",""
"1760","201174","class","Actinobacteria (class)","1","131567","2","","201174","1760","","","","","","",""
"91061","1239","class","Bacilli","1","131567","2","","1239","91061","","","","","","",""
"28216","1224","class","Betaproteobacteria","1","131567","2","","1224","28216","","","","","","",""
"186801","1239","class","Clostridia","1","131567","2","","1239","186801","","","","","","",""
"68337","200930","class","Deferribacteres (class)","1","131567","2","","200930","68337","","","","","","",""
"526524","1239","class","Erysipelotrichi","1","131567","2","","1239","526524","","","","","","",""
"203490","32066","class","Fusobacteria (class)","1","131567","2","","32066","203490","","","","","","",""
"1236","1224","class","Gammaproteobacteria","1","131567","2","","1224","1236","","","","","","",""
"""

def test_taxtable_to_tree():
    h = cStringIO.StringIO(taxtable)
    t = taxtable_to_tree(h)
    assert t.key == "1"
    assert [x.key for x in t.children] == ["131567"]
    assert [x.key for x in t.descendents['2'].children] == ['68336','201174','200930','1239','32066','1224']

def test_lonely_nodes():
    t = Tree(1)(
            Tree(2)(
                Tree(3),
                Tree(4)))
    assert [x.key for x in t.lonelynodes()] == [2]
    t2 = Tree(1)(
             Tree(3)(
                 Tree(4)(
                     Tree(5)(
                         Tree(7)),
                     Tree(6))))
    assert [x.key for x in t2.lonelynodes()] == [3,4,7]
    t2.descendents[5](Tree(9))
    assert 9 in t2.descendents

def test_lonely_company():
    if not(os.path.exists('../testfiles/taxonomy.db')):
        return
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, ncbi.ranks)
    print tax.lineage(49896)
    lonely_tax_ids = [None, 816, 1239]
    for a,t in zip(lonely_tax_ids, lonely_company(tax, lonely_tax_ids)):
        assert a is None or t is None or tax.is_ancestor_of(t, tax.parent_id(a))
        assert t is None or tax.rank(t) == 'species'




