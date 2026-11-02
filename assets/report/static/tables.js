// Sortable, filterable tables for the offline report (no library): click a header to sort, type to filter.
document.querySelectorAll("table.data").forEach(function (table) {
  var tbody = table.tBodies[0], rows = Array.from(tbody.rows);
  var box = document.createElement("input");
  box.type = "search"; box.placeholder = "Filter " + rows.length + " rows"; box.className = "filter";
  var count = document.createElement("span"); count.className = "count";
  table.parentNode.insertBefore(box, table); table.parentNode.insertBefore(count, table);
  box.addEventListener("input", function () {
    var q = box.value.toLowerCase(), shown = 0;
    rows.forEach(function (r) { var hit = r.textContent.toLowerCase().indexOf(q) >= 0; r.hidden = !hit; shown += hit; });
    count.textContent = q ? shown + " of " + rows.length : "";
  });
  Array.from(table.tHead.rows[0].cells).forEach(function (th, i) {
    th.addEventListener("click", function () {
      var dir = th.dataset.dir === "asc" ? "desc" : "asc";
      Array.from(table.tHead.rows[0].cells).forEach(function (c) { delete c.dataset.dir; });
      th.dataset.dir = dir;
      var num = function (s) { var v = parseFloat(s.replace(/,/g, "")); return isNaN(v) ? null : v; };
      rows.sort(function (a, b) {
        var x = a.cells[i].textContent, y = b.cells[i].textContent, nx = num(x), ny = num(y);
        var c = (nx !== null && ny !== null) ? nx - ny : x.localeCompare(y, undefined, {numeric: true});
        return dir === "asc" ? c : -c;
      });
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });
});
