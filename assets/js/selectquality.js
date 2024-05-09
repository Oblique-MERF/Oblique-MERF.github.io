function start() {
  const params = new URL(window.location.href).searchParams;
  dirUrl = '../data/' + params.get('scene');
  const qualityPresets = ['phone', 'low', 'medium', 'high'];
  for (const quality of qualityPresets) {
    console.log(quality);
    const e = document.getElementById(quality);
    e.setAttribute('href', "/viewer?dir=" + dirUrl + '&quality=' + quality);
  }
}

window.onload = start;