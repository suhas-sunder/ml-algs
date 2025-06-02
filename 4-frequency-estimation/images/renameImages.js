// renameImages.js
const fs = require('fs');
const path = require('path');

const folderPath = './';

fs.readdir(folderPath, (err, files) => {
  if (err) return console.error('Error reading directory:', err);

  files.forEach((file) => {
    const match = file.match(/\d.*/); // Find first digit and everything after
    if (match) {
      const newName = match[0];
      const oldPath = path.join(folderPath, file);
      const newPath = path.join(folderPath, newName);

      // Avoid overwriting existing files
      if (fs.existsSync(newPath)) {
        console.warn(`Skipping: ${newName} already exists.`);
        return;
      }

      fs.rename(oldPath, newPath, (err) => {
        if (err) console.error(`Failed to rename ${file}:`, err);
        else console.log(`Renamed: ${file} → ${newName}`);
      });
    } else {
      console.log(`No digits found in: ${file}`);
    }
  });
});
