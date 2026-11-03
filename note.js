(function () {
  var n = document.getElementById('note'), b = document.getElementById('bang'), x = document.getElementById('x');
  function open() { n.classList.add('on'); x.focus(); }
  function close() { n.classList.remove('on'); }
  b.addEventListener('click', open);
  x.addEventListener('click', close);
  addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
})();
