// Eukore's website: the menu on a phone, the charts' tooltips, the statistics page's taxonomy tree (data/taxonomy.json:
// [name, rank, species, with a genome, children | accession], rank an index of its "ranks", -1 the root).
(function () {
  'use strict';

  // --- the menu on a narrow screen ---------------------------------------------------------------------------
  var toggle = document.querySelector('.menu-toggle');
  var nav = document.getElementById('site-menu');
  if (toggle && nav) {
    toggle.addEventListener('click', function () {
      var open = nav.classList.toggle('open');
      toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
    });
  }

  // --- tooltips: any mark with a data-tip ------------------------------------------------------------------
  var tip = document.getElementById('tip');
  function showTip(target, x, y) {
    if (!tip) return;
    tip.textContent = target.getAttribute('data-tip');
    tip.hidden = false;
    var w = tip.offsetWidth, h = tip.offsetHeight;
    var left = Math.min(x + 14, window.innerWidth - w - 8);
    var top = y + 18 + h > window.innerHeight ? y - h - 12 : y + 18;
    tip.style.left = Math.max(8, left) + 'px';
    tip.style.top = Math.max(8, top) + 'px';
  }
  function hideTip() { if (tip) tip.hidden = true; }
  document.addEventListener('mousemove', function (e) {
    var t = e.target.closest && e.target.closest('[data-tip]');
    if (t) showTip(t, e.clientX, e.clientY); else hideTip();
  });
  document.addEventListener('focusin', function (e) {
    var t = e.target.closest && e.target.closest('[data-tip]');
    if (t) { var r = t.getBoundingClientRect(); showTip(t, r.left, r.bottom); }
  });
  document.addEventListener('focusout', hideTip);
  window.addEventListener('scroll', hideTip, { passive: true });

  // --- the taxonomy tree ---------------------------------------------------------------------------------
  var holder = document.getElementById('taxonomy-tree');
  if (!holder) return;
  var PLURAL = { supergroup: 'supergroups', group: 'groups', phylum: 'phyla', 'class': 'classes', order: 'orders',
                 family: 'families', genus: 'genera', species: 'species' };
  var SHOWN = 60;  // children drawn at first; the rest on asking
  var ranks = [];
  var nextId = 0;
  var fmt = new Intl.NumberFormat('en-GB');

  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function kids(node) { return Array.isArray(node[4]) ? node[4] : null; }
  function rankName(node) { return node[1] >= 0 ? ranks[node[1]] : 'root'; }

  function row(node, parent) {
    var li = el('li');
    li.treeNode = node;
    var div = el('div', 'tree-row');
    var children = kids(node);
    if (children && children.length) {
      var b = el('button', 'tree-toggle', '▸');
      var id = 'tree-' + (nextId++);
      b.type = 'button';
      b.setAttribute('aria-expanded', 'false');
      b.setAttribute('aria-controls', id);
      b.setAttribute('aria-label', 'Open ' + node[0]);
      b.addEventListener('click', function () { setOpen(li, b.getAttribute('aria-expanded') !== 'true'); });
      li.listId = id;
      div.appendChild(b);
    } else {
      div.appendChild(el('span'));
    }
    var isSpecies = rankName(node) === 'species';
    var name = el('span', 'tree-name' + (isSpecies ? ' species' : ''), node[0]);
    name.title = node[0];
    div.appendChild(name);
    div.appendChild(el('span', 'tree-rank', rankName(node)));
    var bar = el('span', 'tree-bar');
    var whole = Math.max(1, parent[2]);
    var g = el('span', 'g'), r = el('span', 'r');
    g.style.width = (100 * node[3] / whole) + '%';
    r.style.width = (100 * (node[2] - node[3]) / whole) + '%';
    if (node[3] > 0) bar.appendChild(g);
    if (node[2] - node[3] > 0) bar.appendChild(r);
    var share = Math.round(1000 * node[2] / whole) / 10;
    bar.setAttribute('data-tip', node[0] + ': ' + fmt.format(node[2]) + ' of ' + parent[0] + "'s " +
      fmt.format(parent[2]) + ' species (' + share + '%), ' + fmt.format(node[3]) + ' with a genome or transcriptome');
    div.appendChild(bar);
    var count;
    if (isSpecies) {
      count = el('span', 'tree-count');
      var acc = typeof node[4] === 'string' && node[4] !== '-' ? node[4] : 'rRNA records';
      count.appendChild(el('span', 'tree-acc', acc));
    } else {
      var below = children ? children.length : 0;
      var childRank = below ? rankName(children[0]) : '';
      count = el('span', 'tree-count', fmt.format(node[2]) + ' species' +
        (below && childRank !== 'species' ? ' · ' + fmt.format(below) + ' ' + (below === 1 ? childRank : PLURAL[childRank] || childRank) : ''));
    }
    div.appendChild(count);
    li.appendChild(div);
    return li;
  }

  function fill(ul, parent, limit) {
    var children = kids(parent) || [];
    var n = Math.min(children.length, limit);
    for (var i = ul.childElementCount; i < n; i++) ul.appendChild(row(children[i], parent));
    var more = ul.parentNode && ul.parentNode.querySelector(':scope > .tree-more');
    if (more) more.remove();
    if (children.length > n) {
      var p = el('div', 'tree-more');
      var b = el('button', null, 'Show all ' + fmt.format(children.length) + ' (' + fmt.format(children.length - n) + ' more)');
      b.type = 'button';
      b.addEventListener('click', function () { fill(ul, parent, children.length); });
      p.appendChild(b);
      ul.parentNode.appendChild(p);
    }
  }

  function setOpen(li, open) {
    var b = li.querySelector(':scope > .tree-row > .tree-toggle');
    if (!b) return;
    var ul = li.querySelector(':scope > ul');
    if (open && !ul) {
      ul = el('ul');
      ul.id = li.listId;
      li.appendChild(ul);
      fill(ul, li.treeNode, SHOWN);
    }
    if (ul) ul.hidden = !open;
    var more = li.querySelector(':scope > .tree-more');
    if (more) more.hidden = !open;
    b.setAttribute('aria-expanded', open ? 'true' : 'false');
    b.setAttribute('aria-label', (open ? 'Close ' : 'Open ') + li.treeNode[0]);
  }

  // every taxon with its path, for the search
  var index = [];
  function walk(node, path) {
    var children = kids(node);
    if (node[1] >= 0) index.push({ key: node[0].toLowerCase(), node: node, path: path });
    if (children) for (var i = 0; i < children.length; i++) walk(children[i], path.concat([i]));
  }

  // opens the tree along a path of child indices, then marks and shows the taxon there
  function reveal(rootUl, root, path) {
    var ul = rootUl, node = root, li = null;
    for (var d = 0; d < path.length; d++) {
      var i = path[d];
      if (ul.childElementCount <= i) fill(ul, node, i + 1);
      li = ul.children[i];
      node = li.treeNode;
      if (d < path.length - 1) { setOpen(li, true); ul = li.querySelector(':scope > ul'); }
    }
    if (!li) return;
    var old = holder.querySelector('.tree-row.found');
    if (old) old.classList.remove('found');
    var r = li.querySelector(':scope > .tree-row');
    r.classList.add('found');
    r.scrollIntoView({ block: 'center', behavior: 'smooth' });
  }

  fetch(holder.getAttribute('data-src'))
    .then(function (resp) { if (!resp.ok) throw new Error(resp.status); return resp.json(); })
    .then(function (data) {
      ranks = data.ranks;
      var root = data.tree;
      holder.textContent = '';
      var ul = el('ul');
      holder.appendChild(ul);
      fill(ul, root, Infinity);
      // the supergroups open at first: their groups in view
      for (var i = 0; i < ul.children.length && i < 1; i++) setOpen(ul.children[i], true);
      walk(root, []);
      var search = document.getElementById('tree-search');
      var results = document.getElementById('tree-results');
      if (!search || !results) return;
      search.addEventListener('input', function () {
        var q = search.value.trim().toLowerCase();
        results.textContent = '';
        if (q.length < 2) { results.hidden = true; return; }
        var hits = [];
        for (var j = 0; j < index.length && hits.length < 40; j++) {
          if (index[j].key.indexOf(q) >= 0) hits.push(index[j]);
        }
        hits.sort(function (a, b) { return (a.key.indexOf(q) !== 0) - (b.key.indexOf(q) !== 0) || a.node[1] - b.node[1]; });
        if (!hits.length) results.appendChild(el('p', null, 'No taxon of this release is named so.'));
        hits.forEach(function (h) {
          var b = el('button', null, h.node[0] + ' ');
          b.type = 'button';
          b.appendChild(el('span', null, rankName(h.node) + ', ' + fmt.format(h.node[2]) + ' species'));
          b.addEventListener('click', function () { results.hidden = true; reveal(ul, root, h.path); });
          results.appendChild(b);
        });
        results.hidden = false;
      });
    })
    .catch(function () {
      holder.insertAdjacentHTML('beforeend', '<p class="notice">The taxonomy could not be loaded (data/taxonomy.json).</p>');
    });
})();
